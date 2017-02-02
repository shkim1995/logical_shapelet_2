 
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

using namespace std; 


class Type{
public:
	int type;
	int count;

	Type(){
		type = 0;
		count = 0;
	}

	Type(int t, int c){
		type = t;
		count = c;
	}

	void print(){
		cout<<type<<" : "<<count<<endl;
	}
};


//match data type and number of data
class TypeCount{
public:
	list<Type> types;
	int typeNum;
	int totalCount;

	TypeCount(){
		typeNum = 0;
		totalCount = 0;
	}

	TypeCount(vector<Data>* dataVector){
		totalCount = 0;
		for(int i=0; i<dataVector->size(); i++){
			insert((*dataVector)[i].type);
		}
	}

	int count(int type){
		list<Type>::iterator it = types.begin();
		while(it!=types.end()){
			if(it->type==type)
				break;
			it++;
		}
		if(it!=types.end())
			return it->count;
		return -1;
	}


	void insert(int newType){
		typeNum = 0;
		list<Type>::iterator it = types.begin();
		while(true){
			if(it==types.end()||it->type>newType){
				Type t(newType, 1);
				types.insert(it, t);
				typeNum++;
				totalCount++;
				break;
			}
			else if(it->type==newType){
				it->count++;
				totalCount++;
				break;
			}
			it++;
		}
	}

	void add(int type){
		totalCount++;

		list<Type>::iterator it = types.begin();
		while(it!=types.end()){
			if(it->type==type)
				break;
			it++;
		}
		if(it!=types.end())
			it->count++;
	}

	void subtract(int type){
		totalCount--;

		list<Type>::iterator it = types.begin();
		while(it!=types.end()){
			if(it->type==type)
				break;
			it++;
		}
		if(it!=types.end()){
			it->count--;
		}
	}

	void print(){

		list<Type>::iterator it = types.begin();

		while(it!=types.end()){
			it->print();
			it++;
		}
	}

	TypeCount defaultCount(){
		TypeCount newTypeCount;
		list<Type>::iterator it = types.begin();
		while(it!=types.end()){
			Type temp(it->type, 0);
			newTypeCount.types.push_back(temp);
			newTypeCount.typeNum++;
			it++;
		}
		return newTypeCount;
	}
};
 
