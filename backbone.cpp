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
	for(auto it(kmerCount.begin());it!=kmerCount.end();++it){
		if(hasMom.count(it->first)==0){
			result.push_back(it->first);
			kmerCount.erase(it->first);
		}
	}
	return result;
}



vector<string> noDaughter(vector<pair<string,string>>& pairList, unordered_map<string,uint>& kmerCount){
	unordered_set<string> hasMom;
	vector<string> result;
	for(uint i(0);i<pairList.size();++i){
		if(pairList[i].first!=""){
			hasMom.insert(pairList[i].first);
		}
	}
	for(uint i(0);i<pairList.size();++i){
		if(pairList[i].first!=""){
			if(hasMom.count(pairList[i].second)==0){
				pairList[i]={"",""};
			}
		}
	}
	for(auto it(kmerCount.begin());it!=kmerCount.end();++it){
		if(hasMom.count(it->first)==0){
			result.push_back(it->first);
			kmerCount.erase(it->first);
		}
	}
	return result;
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
	//READ LOADING
	while(not in.eof()){
		getline(in,header);
		getline(in,read);
		readsVector.push_back(read);
	}
	//KMER COUNTING
	unordered_map<string,uint> kmerCount;
	for(uint i(0);i<readsVector.size();++i){
		read=readsVector[i];
		for(uint ii(0);ii<read.size();++ii){
			kmerCount[read.substr(ii,k)]++;
		}
	}
	for(auto it(kmerCount.begin()); it != kmerCount.end(); ++it){
		if(it->second<solid){
			kmerCount.erase(it->first);
		}
	}
	//SKETCH CREATION
	string kmer;
	vector<vector<string>> sketchVector(readsVector.size());
	for(uint i(0);i<readsVector.size();++i){
		read=readsVector[i];
		for(uint ii(0);ii<read.size();++ii){
			kmer=read.substr(ii,k);
			if(kmerCount[kmer]>=solid){
				sketchVector[ii].push_back(kmer);
			}
		}
	}
	vector<string> sketch;
	vector<pair<string,string>> pairList;
	//PAIR LIST CREATION
	for(uint i(0);i<sketchVector.size();++i){
		sketch=sketchVector[i];
		for(uint ii(0);ii<sketch.size()-1;++ii){
			pairList.push_back({sketch[i],sketch[i+1]});
		}
	}
	vector<string> resultBegin, resultEnd,nomom,nodaughter,result;
	while(not kmerCount.empty()){
		nomom=noMom(pairList,kmerCount);
		nodaughter=noDaughter(pairList,kmerCount);
		if(nomom.size()==1){
			resultBegin.push_back(nomom[0]);
		}
		if(nodaughter.size()==1){
			resultEnd.push_back(nomom[0]);
		}
	}
	reverse(resultEnd.begin(), resultEnd.end());
	resultBegin.insert(resultBegin.end(), resultEnd.begin(), resultEnd.end());

	for(uint i(0);i<resultBegin.size();++i){
		cout<<resultBegin[i]<<" ";
	}
	cout<<endl;

	return 0;
}

