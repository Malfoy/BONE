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



vector<string> noMom(vector<pair<string,string>>& pairList, unordered_map<string,uint>& kmerCount){
	unordered_set<string> hasMom;
	vector<string> result;
	for(uint i(0);i<pairList.size();++i){
		if(pairList[i].first!=""){
			hasMom.insert(pairList[i].second);
		}
	}
	for(uint i(0);i<pairList.size();++i){
		if(pairList[i].first!=""){
			if(hasMom.count(pairList[i].first)==0){
				pairList[i]={"",""};
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



vector<string> noDaughter(vector<pair<string,string>>& pairList, unordered_map<string,uint>& kmerCount){
	unordered_set<string> hasDaughter;
	vector<string> result;
	//~ cout<<"lol1"<<endl;
	for(uint i(0);i<pairList.size();++i){
		if(pairList[i].first!=""){
			hasDaughter.insert(pairList[i].first);
		}
	}
	//~ cout<<"lol2"<<endl;
	for(uint i(0);i<pairList.size();++i){
		if(pairList[i].first!=""){
			if(hasDaughter.count(pairList[i].second)==0){
				pairList[i]={"",""};
			}
		}
	}
	//~ cout<<"lol3"<<endl;
	vector<string> toErase;
	for(auto it(kmerCount.begin());it!=kmerCount.end();++it){
		if(hasDaughter.count(it->first)==0){
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



void eraseKmer(string kmer, vector<pair<string,string>>& pairList, unordered_map<string,uint>& kmerCount){
	for(uint i(0);i<pairList.size();++i){
		if(pairList[i].first!=""){
			if(kmer==pairList[i].second){
				pairList[i]={"",""};
			}
		}
	}

	kmerCount.erase(kmer);
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
	vector<pair<string,string>> pairList;
	//PAIR LIST CREATION
	for(uint i(0);i<sketchVector.size();++i){
		sketch=sketchVector[i];
		for(uint ii(0);ii+1<sketch.size();++ii){
			pairList.push_back({sketch[ii],sketch[ii+1]});
		}
	}
	cout<<4<<endl;
	vector<string> resultBegin, resultEnd,nomom,nodaughter,result;
	uint indiceSupress(0);
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
			unordered_set<string> kmersToSupress;
			sketch=sketchVector[indiceSupress%sketchVector.size()];
			indiceSupress++;
			for(uint ii(0);ii<sketch.size();++ii){
				if(kmerCount.count(sketch[ii])!=0){
					eraseKmer(sketch[ii],pairList,kmerCount);
					ii=sketch.size();
				}
			}
			cout<<"left"<<kmerCount.size()<<endl;
		}else{
			indiceSupress=0;
		}
	}
	cout<<5<<endl;
	reverse(resultEnd.begin(), resultEnd.end());
	resultBegin.insert(resultBegin.end(), resultEnd.begin(), resultEnd.end());
	//~ cout<<"Result"<<endl;
	unordered_map<string, uint> kmersBackbone;
	for(uint i(0); i < resultBegin.size(); ++i){
		cout<< resultBegin[i] << " ";
		kmersBackbone.insert({resultBegin[i], i});
	}
	cout<<endl;

	//~ string kmer;
	for(uint i(0); i < readsVector.size(); ++i){
		read = readsVector[i];
		for(uint ii(0); ii + k <= read.size(); ++ii){
			kmer = read.substr(ii,k);
			if (kmersBackbone.count(kmer)){
				cout << kmersBackbone[kmer] << " ";
			}
		}
		cout << endl;
	}

	return 0;
}


