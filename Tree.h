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


#ifndef INFO
#define INFO
#include "Info.h"
#endif

using namespace std;


class Node{

public:

	ShapeletInfo shapeletInfo;
	int type; // type at the bottom node, -1 for parents

	Node* left;
	Node* right;

	Node(){

		type = -1;
		left = NULL;
		right = NULL;

	}

	void printNodes(){
		if(type==-1)
			cout<<shapeletInfo.best_t<<endl;
		else
			cout<<type<<endl;
		//printVector(shapelet);
		if(left!=NULL)
			left->printNodes();

		if(right!=NULL)
			right->printNodes();
	}

private:

	void printVector(vector<double> S){
	for(int i=0; i<S.size(); i++)
		cout<<S[i]<<" ";
	cout<<endl;
	}

};
