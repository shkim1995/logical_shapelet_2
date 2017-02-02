#ifndef HEADERS

#define HEADERS
#include <iostream>
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

using namespace std;

//individual item for Histogram
//including type and value
class HistData{
public:

	int type;
	double value;	

	//CONSTRUCTORS//
	HistData(){
		type = 0;
		value = 0;
	}

	HistData(int t, double v){
		type = t;
		value = v;
	}
	/////////////////


	void print(){
		if(value<0.0001)	
			cout<<type<<":"<<0<<" | ";
		else
			cout<<type<<":"<<value<<" | ";
	}

};


//Histogram - Automatically Sorted by Value
class Histogram{
public:
	list<HistData> data;

	void insert(int type, double value){
		HistData newData(type, value);
		if(data.size()==0){
			data.push_back(newData);
			return;
		}

		list<HistData>::iterator it = data.begin();
		while(it!=data.end()){
			if(it->value>value)
				break;
			it++;
		}
		data.insert(it, newData);
	}

	list<HistData>::iterator begin(){
		return this->data.begin();
	}

	list<HistData>::iterator end(){
		return this->data.end();
	}

	int size(){
		return data.size();
	}


	void print(){
		if(this==NULL)
			return;

		list<HistData>::iterator it = data.begin();
		while(it!=data.end()){
			it->print();
			it++;
		}
		cout<<endl;
	}

};