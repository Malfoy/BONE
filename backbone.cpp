#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>



using namespace std;



char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}



string revComp(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = revCompChar(s[i]);
	}
	return rc;
}



string getCanonical(const string& str){
	return (min(str,revComp(str)));
}



vector<string> withoutAdj(unordered_map<string,vector<string>>& prev, unordered_map<string,vector<string>>& next, unordered_map<string,uint>& kmerCount, vector<string>& candidate){
	vector<string> nextCandidate,res,toSuppr;
	string kmer,son;
	if(candidate.size()==0){
		for(auto it(kmerCount.begin());it!=kmerCount.end();++it){
			kmer=(it->first);
			if(prev.count(kmer)==0 or prev[kmer].size()==0){
				//~ cout<<"noprevfound"<<endl;
 				toSuppr.push_back(kmer);
				res.push_back(kmer);
				for(uint j(0);j<next[kmer].size();++j){
					son=next[kmer][j];
					nextCandidate.push_back(son);
					auto newEnd = std::remove(prev[son].begin(), prev[son].end(), kmer);
					prev[son].erase(newEnd, prev[son].end());
				}
				next.erase(kmer);
				prev.erase(kmer);
			}
		}
		for(uint i(0);i<toSuppr.size();++i){
			kmerCount.erase(toSuppr[i]);
		}
	}else{
		for(uint i(0);i<candidate.size();++i){
			kmer=candidate[i];
			if(kmerCount.count(kmer)!=0){
				if(prev.count(kmer)==0 or prev[kmer].size()==0){
					//~ cout<<"noprevfound"<<endl;
					kmerCount.erase(kmer);
					res.push_back(kmer);
					for(uint j(0);j<next[kmer].size();++j){
						son=next[kmer][j];
						//~ if(kmerCount.count(son)!=0){
							nextCandidate.push_back(son);
						//~ }
						auto newEnd = std::remove(prev[son].begin(), prev[son].end(), kmer);
						prev[son].erase(newEnd, prev[son].end());
					}
					next.erase(kmer);
					prev.erase(kmer);
				}
			}else{
				//~ cout<<"wow"<<endl;
			}
		}
	}
	sort( nextCandidate.begin(), nextCandidate.end() );
	nextCandidate.erase( unique( nextCandidate.begin(), nextCandidate.end() ), nextCandidate.end() );
	candidate=nextCandidate;
	return res;
}



void eraseKmer(unordered_map<string,vector<string>>& prev, unordered_map<string,vector<string>>& next, unordered_map<string,uint>& kmerCount,string kmer){
	kmerCount.erase(kmer);
	string son;
	for(uint j(0);j<next[kmer].size();++j){
		son=next[kmer][j];
		auto newEnd = std::remove(prev[son].begin(), prev[son].end(), kmer);
		prev[son].erase(newEnd, prev[son].end());
	}
	next.erase(kmer);
	for(uint j(0);j<prev[kmer].size();++j){
		son=prev[kmer][j];
		auto newEnd = std::remove(next[son].begin(), next[son].end(), kmer);
		next[son].erase(newEnd, next[son].end());
	}
	prev.erase(kmer);
}



unordered_set<string> getFirstKmers(uint depth, unordered_map<string,uint>& kmerCount, vector<vector<string>>& sketchVector){
	unordered_set<string> res;
	for(uint i(0);i<sketchVector.size();++i){
		uint n(0);
		for(uint j(0);j<sketchVector[i].size();++j){
			if(kmerCount.count(sketchVector[i][j])==1){
				res.insert(sketchVector[i][j]);
				if(++n>=depth){
					break;
				}
			}
		}
	}
	//~ cout<<res.size()<<endl;
	return res;
}



unordered_set<string> getLastKmers(uint depth, unordered_map<string,uint>& kmerCount, vector<vector<string>>& sketchVector){
	unordered_set<string> res;
	for(uint i(0);i<sketchVector.size();++i){
		uint n(0);
		for(int j(sketchVector[i].size()-1);j>=0;--j){
			if(kmerCount.count(sketchVector[i][j])==1){
				res.insert(sketchVector[i][j]);
				if(++n>=depth){
					break;
				}
			}
		}
	}
	//~ cout<<res.size()<<endl;
	return res;
}




//TODO RC ?
int main(int argc, char ** argv){
	if(argc<4){
		cout<<"[read file (fasta oneline)] [kmer size] [solidity threshold value]"<<endl;
		exit(0);
	}
	ifstream in(argv[1]);
	uint k(stoi(argv[2]));
	uint solid(stoi(argv[3]));
	string query,read,header;
	getline(in,header);
	getline(in,query);
	vector<string> readsVector;
	readsVector.push_back(query);
	//READ LOADING
	while(not in.eof()){
		getline(in,header);
		getline(in,read);
		if(read.size()>=k){
			readsVector.push_back(read);
		}
		read="";
	}

	//KMER COUNTING
	unordered_map<string,uint> kmerCount,kmerCountLocal;
	for(uint i(0);i<readsVector.size();++i){
		read=readsVector[i];
		kmerCountLocal={};
		for(uint ii(0);ii+k<=read.size();++ii){
			kmerCountLocal[read.substr(ii,k)]++;
		}
		for(auto it(kmerCountLocal.begin()); it != kmerCountLocal.end(); ++it){
			if(it->second==1){
				kmerCount[it->first]++;
			}
		}
	}
	vector<string> toErase;
	for(auto it(kmerCount.begin()); it != kmerCount.end(); ++it){
		if(it->second<solid){
			toErase.push_back(it->first);
		}
	}
	for(uint i(0);i<toErase.size();++i){
		kmerCount.erase(toErase[i]);
	}

	//SKETCH CREATION
	string kmer;
	vector<vector<string>> sketchVector(readsVector.size());
	for(uint i(0);i<readsVector.size();++i){
		read=readsVector[i];
		for(uint ii(0);ii+k<=read.size();++ii){
			kmer=read.substr(ii,k);
			if(kmerCount.count(kmer)==1){
				sketchVector[i].push_back(kmer);
			}
		}
	}

	vector<string> sketch;
	unordered_map<string,vector<string>> next,prev;
	//ADJ CREATION
	for(uint i(0);i<sketchVector.size();++i){
		sketch=sketchVector[i];
		for(uint ii(0);ii+1<sketch.size();++ii){
			next[sketch[ii]].push_back(sketch[ii+1]);
			prev[sketch[ii+1]].push_back(sketch[ii]);
		}
	}
	vector<string> resultBegin, resultEnd, noPrev, noNext, result,candidateNext,candidatePrev;
	unordered_set<string> kmersBegin,kmersEnd;
	//MAIN LOOP
	while(not kmerCount.empty()){
		//~ cout<<"go"<<endl;
		uint count(kmerCount.size());
		noPrev=withoutAdj(prev, next,kmerCount,candidatePrev);
		noNext=withoutAdj(next, prev,kmerCount,candidateNext);
		if(noPrev.size()==1){
			resultBegin.push_back(noPrev[0]);
		}
		if(noNext.size()==1){
			resultEnd.push_back(noNext[0]);
		}
		//IF PARADOX
		if(kmerCount.size()==count){
			//~ cout<<"paradox"<<endl;cin.get();
			bool found(false);
			uint depth(1);
			while(not found){
				if(kmerCount.size()<2*depth){
					kmerCount={};
					break;
				}
				kmersBegin=getFirstKmers(depth,kmerCount,sketchVector);
				kmersEnd=getLastKmers(depth,kmerCount,sketchVector);
				for(auto it(kmersBegin.begin()); it != kmersBegin.end() and not found; ++it){
					if(kmersEnd.count(*it)==1){
						eraseKmer(prev,next,kmerCount,*it);
						found=true;
						//~ break;
					}
				}
				candidatePrev=candidateNext={};
				depth*=2;
			}

		}
	}
	cout<<resultBegin.size()<<" "<<resultEnd.size()<<endl;
	ofstream out("outBone.txt");
	reverse(resultEnd.begin(), resultEnd.end());
	resultBegin.insert(resultBegin.end(), resultEnd.begin(), resultEnd.end());
	unordered_map<string, uint> kmersBackbone;
	for(uint i(0); i < resultBegin.size(); ++i){
		kmersBackbone.insert({resultBegin[i], i});
	}
	for(uint i(0); i < readsVector.size(); ++i){
		read = readsVector[i];
		int last(-1);
		for(uint ii(0); ii + k <= read.size(); ++ii){
			kmer = read.substr(ii,k);
			if (kmersBackbone.count(kmer)){
				uint indice(kmersBackbone[kmer]);
				if((int)indice<last){
					kmersBackbone.erase(kmer);
				}else{
					last=indice;
				}
			}
		}
	}
	uint countNew(0);
	for(uint i(0); i < resultBegin.size(); ++i){
		if(kmersBackbone.count(resultBegin[i])==1){
			kmersBackbone[resultBegin[i]]=countNew++;
			out<< resultBegin[i] << " ";
		}
	}
	out<<endl;
	for(uint i(0); i < readsVector.size(); ++i){
		read = readsVector[i];
		for(uint ii(0); ii + k <= read.size(); ++ii){
			kmer = read.substr(ii,k);
			if (kmersBackbone.count(kmer)){
				out << kmersBackbone[kmer] << " ";
			}
		}
		out << endl;
	}



	return 0;
}


