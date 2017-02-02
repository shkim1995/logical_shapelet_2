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
 
using namespace std;


class Data{
	public:
		int type;
		vector<double> data;
		int type2;
 
		void print(){
			cout<<type<<" | ";
			for(int i=0; i<data.size(); i++)
				cout<<" "<<data[i];
			cout<<endl;
			cout<<endl;
		}
 
		int size(){
			return data.size();
		}
};