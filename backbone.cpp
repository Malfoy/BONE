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



vector<string> noMom(vector<vector<pair<string,string>>>& pairList, unordered_map<string,uint>& kmerCount){
	unordered_set<string> hasMom;
	vector<string> result;
	for(uint j(0);j<pairList.size();++j){
		for(uint i(0);i<pairList[j].size();++i){
			if(pairList[j][i].first!=""){
				hasMom.insert(pairList[j][i].second);
			}
		}
	}
	for(uint j(0);j<pairList.size();++j){
		for(uint i(0);i<pairList[j].size();++i){
			if(pairList[j][i].first!=""){
				if(hasMom.count(pairList[j][i].first)==0){
					pairList[j][i]={"",""};
				}
			}
		}
	}
	vector<string> toErase;
	for(auto it(kmerCount.begin());it!=kmerCount.end();++it){
		if(hasMom.count(it->first)==0){
			//~ cout<<it->first<<endl;
			result.push_back(it->first);
			toErase.push_back(it->first);
			//~ kmerCount.erase(it->first);
		}
	}
	for(uint i(0);i<toErase.size();++i){
		kmerCount.erase(toErase[i]);
	}
	return result;
}



vector<string> noDaughter(vector<vector<pair<string,string>>>& pairList, unordered_map<string,uint>& kmerCount){
	unordered_set<string> hasDaughter;
	vector<string> result;
	for(uint j(0);j<pairList.size();++j){
		for(uint i(0);i<pairList[j].size();++i){
			if(pairList[j][i].first!=""){
				hasDaughter.insert(pairList[j][i].first);
			}
		}
	}
	for(uint j(0);j<pairList.size();++j){
		for(uint i(0);i<pairList[j].size();++i){
			if(pairList[j][i].first!=""){
				if(hasDaughter.count(pairList[j][i].second)==0){
					pairList[j][i]={"",""};
				}
			}
		}
	}
	vector<string> toErase;
	for(auto it(kmerCount.begin());it!=kmerCount.end();++it){
		if(hasDaughter.count(it->first)==0){
			result.push_back(it->first);
			toErase.push_back(it->first);
		}
	}
	for(uint i(0);i<toErase.size();++i){
		kmerCount.erase(toErase[i]);
	}
	return result;
}



void eraseKmer(string kmer, vector<vector<pair<string,string>>>& pairList, unordered_map<string,uint>& kmerCount){
	for(uint j(0);j<pairList.size();++j){
		for(uint i(0);i<pairList[j].size();++i){
			if(pairList[j][i].first!=""){
				if(kmer==pairList[j][i].second or kmer==pairList[j][i].first){
					pairList[j][i]={"",""};
				}
			}
		}
	}
	kmerCount.erase(kmer);
}



unordered_set<string> getFirstKmers(uint depth, vector<vector<pair<string,string>>>& pairList){
	unordered_set<string> res;
	for(uint i(0);i<pairList.size();++i){
		uint n(0);
		for(uint j(0);j<pairList[i].size() and n<=depth;++j){
			if(pairList[i][j].first!=""){
				res.insert(pairList[i][j].first);
				++n;
			}
		}
	}
	return res;
}



unordered_set<string> getLastKmers(uint depth, vector<vector<pair<string,string>>>& pairList){
	unordered_set<string> res;
	for(uint i(0);i<pairList.size();++i){
		uint n(0);
		for(int j(pairList[i].size()-1);j>=0 and n<=depth;--j){
			if(pairList[i][j].first!=""){
				res.insert(pairList[i][j].first);
				++n;
			}
		}
	}
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
	cout<<1<<endl;
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
	cout<<3<<endl;
	vector<string> sketch;
	vector< vector<pair<string,string>> > pairList(sketchVector.size());
	//PAIR LIST CREATION
	for(uint i(0);i<sketchVector.size();++i){
		sketch=sketchVector[i];
		for(uint ii(0);ii+1<sketch.size();++ii){
			pairList[i].push_back({sketch[ii],sketch[ii+1]});
		}
	}
	cout<<4<<endl;
	vector<string> resultBegin, resultEnd,nomom,nodaughter,result;
	unordered_set<string> kmersBegin,kmersEnd;
	while(not kmerCount.empty()){
		uint count(kmerCount.size());
		cout<<kmerCount.size()<<endl;
		nomom=noMom(pairList,kmerCount);
		nodaughter=noDaughter(pairList,kmerCount);
		if(nomom.size()==1){
			resultBegin.push_back(nomom[0]);
		}
		if(nodaughter.size()==1){
			resultEnd.push_back(nodaughter[0]);
		}
		if(kmerCount.size()==count){
			cout<<"go"<<endl;
			bool found(false);
			uint depth(0);
			while(not found){
				if(kmerCount.size()<2*depth){
					kmerCount={};
				}
				cout<<"cacao"<<endl;
				kmersBegin=getFirstKmers(depth,pairList);
				cout<<"gogogo"<<endl;
				kmersEnd=getLastKmers(depth,pairList);
				cout<<"chocolat"<<endl;
				for(auto it(kmersBegin.begin()); it != kmersBegin.end(); ++it){
					if(kmersEnd.count(*it)==1){
						eraseKmer(*it,pairList,kmerCount);
						found=true;
					}
				}
				++depth;

			}

		}
	}
	cout<<5<<endl;
	ofstream out("outBone.txt");
	reverse(resultEnd.begin(), resultEnd.end());
	resultBegin.insert(resultBegin.end(), resultEnd.begin(), resultEnd.end());
	//~ cout<<"Result"<<endl;
	unordered_map<string, uint> kmersBackbone;
	for(uint i(0); i < resultBegin.size(); ++i){
		out<< resultBegin[i] << " ";
		kmersBackbone.insert({resultBegin[i], i});
	}
	out<<endl;

	//~ string kmer;
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


