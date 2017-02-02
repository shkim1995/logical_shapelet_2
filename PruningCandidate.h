#ifndef HEADERS
 
#define HEADERS
 
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <sstream>
#include <limits>
#include <cstddef>
#include <cstring>
#include "time.h"
 
#endif
 
#ifndef DATA
#define DATA
#include "Data.h"
#endif


#ifndef INFO
#define INFO
#include "Info.h"
#endif


#ifndef HIST
#define HIST
#include "Histogram.h"
#endif

using namespace std;


class PruningCandidateItem{
public:
	Histogram hist;
	int k; //data Number
	int u; //starting point
	int l; // length

	PruningCandidateItem(){

	}

	PruningCandidateItem(Histogram h, int x, int y, int z){
		hist = h;
		k = x;
		u = y;
		l = z;
	}

	void print(){
		cout<<k<<" "<<u<<" "<<l<<" | ";
	}
};



class PruningCandidate{
public:
	list<PruningCandidateItem> candidateList;

	void insert(Histogram* hist, int k, int u, int l){
		PruningCandidateItem temp(*hist, k, u, l);
		candidateList.push_back(temp);
		
		// cout<<"DEBUG1"<<endl;
		// print();

		if(candidateList.size()>CANLEN){
			list<PruningCandidateItem>::iterator it = candidateList.begin();
			candidateList.erase(it);
		}

	}

	void clear(){
		while(true){
			if(candidateList.size()==0)
				break;
			candidateList.erase(candidateList.begin());
		}
	}

	list<PruningCandidateItem>::iterator begin(){
		return candidateList.begin();
	}

	list<PruningCandidateItem>::iterator end(){
		return candidateList.end();
	}

	void print(){
		cout<<"length "<<candidateList.size()<<" : ";
		list<PruningCandidateItem>::iterator it = candidateList.begin();
		while(it!=candidateList.end()){
			it->print();
			it++;
		}
		cout<<endl;

	}
		
};
 