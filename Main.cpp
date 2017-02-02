#include "Stats.h"
#include "Tree.h"
#include "TypeCount.h"
#include "PruningCandidate.h"


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
 
 
//////////////////////////////////////////////////////////////////////////
/////////////////////////////MAIN FUNCTIONS///////////////////////////////
//////////////////////////////////////////////////////////////////////////
 
 
//FOR DEBUGGINH
void printVector(vector<double> S){
	for(int i=0; i<S.size(); i++)
		cout<<S[i]<<" ";
	cout<<endl;
}
 
 
void printShapelet(ShapeletInfo info){
	cout<<"Shapelet : "<<endl;
	printVector(info.best_S);
	cout<<"Distance : "<<info.best_t<<endl;
 
}
 
 
//normalize dataVector[from] ~ dataVector[to-1]
vector<double> zNorm(vector<double>* dataVector, int from, int to){
 
    vector<double> normalizedVector;
 
	double sum = 0;
	int num = 0;
	double avg = 0;
	double var = 0;
 
	for(int i=from; i<to; i++){
		sum = sum + (*dataVector)[i];
		num++;
	}
 
	if(num!=0)	
		avg = sum/num;
	
	for(int i=from; i<to; i++){
		double x = (*dataVector)[i];
		var = var + (x-avg)*(x-avg);
	}
	var = var/(num);
	var = sqrt(var);
 
	for(int i=from; i<to; i++){
		normalizedVector.push_back(((*dataVector)[i]-avg)/var);
	}
 
	return normalizedVector;
	
}
 
//brute force
double sdist(vector<double>* x, vector<double>* y){ // |x|<=|y|
 
	double minSum = 100000;
	double sum;
 
	vector<double> x_norm = zNorm(x, 0, x->size());
	int l = x->size();
 
	for(int i=0; i<y->size()-l+1; i++){
		sum = 0;
		vector<double> y_norm = zNorm(y, i, i+l);
		for(int k=0; k<l; k++){
			sum = sum+(y_norm[k]-x_norm[k])*(y_norm[k]-x_norm[k]);
		}
		//cout<<sqrt(sum/l)<<" ";
		if(minSum>sum)
			minSum = sum;
	}
	//cout<<endl;
	return sqrt(minSum/l);
}
 
double sdist_new(vector<Data>* dataVector, Stats* stats, int i, int j, int u, int l){ 
	// x = Di, y = Dj
	// |x| <= |y|
 
	// cout<<stats->S2[i][0]<<" "<<stats->S2[i][1]<<" "<<stats->S2[i][2]<<" "<<stats->S2[i][3]<<" "<<stats->S2[i][4]<<" ";
	// cout<<stats->S2[i][5]<<" "<<stats->S2[i][6]<<" "<<stats->S2[i][7]<<" "<<stats->S2[i][8]<<" "<<endl;
	// cout<<endl;
 
	// cout<<stats->S2[j][0]<<" "<<stats->S2[j][1]<<" "<<stats->S2[j][2]<<" "<<stats->S2[j][3]<<" "<<stats->S2[j][4]<<" ";
	// cout<<stats->S2[j][5]<<" "<<stats->S2[j][6]<<" "<<stats->S2[j][7]<<" "<<stats->S2[j][8]<<" "<<endl;
	// cout<<endl;
 
	int x = ((*dataVector)[i]).data.size();
	int y = ((*dataVector)[j]).data.size();
 
	double avg_x = (stats->S[i][u+l]-stats->S[i][u])/l;
	double avg_y;
	// cout<<stats->S2[i][u+l]-stats->S2[i][u]<<endl;
	// cout<<endl;
 
	double stdev_x = (stats->S2[i][u+l]-stats->S2[i][u]+0.0)/l - avg_x*avg_x;
	stdev_x = sqrt(stdev_x);
	double stdev_y;
 
	// cout<<"avg_x :"<<avg_x<<endl;
	// cout<<"stdev_x : "<<sqrt(stdev_x)<<endl;
	// cout<<endl;
 
	double multiple;
	double Cs;
	double dist;
	double mindist = 100000;
 
	for(int v=0; v<=y-l; v++){
		avg_y = (stats->S[j][v+l]-stats->S[j][v])/l;
		stdev_y= (stats->S2[j][v+l]-stats->S2[j][v])/l - avg_y*avg_y;
		stdev_y = sqrt(stdev_y);
		// cout<<"avg_y :"<<avg_y<<endl;
		// cout<<"stdev_y : "<<sqrt(stdev_y)<<endl;
 
		multiple = (stats->M[i][j][u+l][v+l]-stats->M[i][j][u][v]);
 
		Cs = (multiple-l*avg_x*avg_y)/(l*stdev_x*stdev_y);
		if(abs(1-Cs)<0.001)
			dist = 0;
		else
			dist = sqrt(2*(1-Cs));
		//cout<<dist<<endl;
 
		//cout<<endl;
 
		if(dist<mindist)
			mindist = dist;
	}
	//cout<<endl;
 
	return mindist;
}

double sdist_new(vector<Data>* dataVector, Stats* stats, int i, int j, int u, int v, int l){ 
	// x = Di, y = Dj
	// |x| <= |y|
 
	int x = ((*dataVector)[i]).data.size();
	int y = ((*dataVector)[j]).data.size();
 
	double avg_x = (stats->S[i][u+l]-stats->S[i][u])/l;
	double avg_y = (stats->S[j][v+l]-stats->S[j][v])/l;
	// cout<<stats->S2[i][u+l]-stats->S2[i][u]<<endl;
	// cout<<endl;
 
	double stdev_x = (stats->S2[i][u+l]-stats->S2[i][u]+0.0)/l - avg_x*avg_x;
	stdev_x = sqrt(stdev_x);
	double stdev_y = (stats->S2[j][v+l]-stats->S2[j][v]+0.0)/l - avg_y*avg_y;
	stdev_y = sqrt(stdev_y);
 
	// cout<<"avg_x :"<<avg_x<<endl;
	// cout<<"stdev_x : "<<sqrt(stdev_x)<<endl;
	// cout<<endl;
 
	double multiple  = (stats->M[i][j][u+l][v+l]-stats->M[i][j][u][v]);
	double Cs = (multiple-l*avg_x*avg_y)/(l*stdev_x*stdev_y);
	double dist;

	if(abs(1-Cs)<0.001)
		dist = 0;
	else
		dist = sqrt(2*(1-Cs));
 
	return dist;
}
 
 
//Entropy before division
// a, b : number of each elements 
double Entropy(int  a, int b){
	if(a==0) return 0;
	if(b==0) return 0;
	double pa = (a+0.0)/(a+b+0.0);
	double pb = (b+0.0)/(a+b+0.0);
	return -(pa*log(pa)+pb*log(pb));
}


double Entropy(TypeCount* typeCount){
	list<Type>::iterator it = typeCount->types.begin();
	double entropy = 0;
	while(it!=typeCount->types.end()){
		if(it->count==0){
			it++;
			continue;
		}
		double p = (it->count+0.0)/(typeCount->totalCount+0.0);
		//cout<<p<<endl;
		entropy -= p*log(p);
		it++;
	}
	return entropy;
}
 
//Entropy after division
double Entropy(int a1, int b1, int a2, int b2){
	
	return ((a1+b1)*Entropy(a1, b1)+(a2+b2)*Entropy(a2, b2))/(a1+a2+b1+b2);
}
 

double Entropy(TypeCount* tc1, TypeCount* tc2){
 	int count1 = tc1->totalCount;
 	int count2 = tc2->totalCount;
 	// cout<<"1 : "<<count1<<", 2: "<<count2<<endl;
 	// cout<<"1 : "<<Entropy(tc1)<<", 2: "<<Entropy(tc2)<<endl;
 	return (count1*Entropy(tc1)+count2*Entropy(tc2))/(count1+count2);
 } 
 


double Gap(double distSum_left, double distSum_right, int num_left, int num_right){
	if(num_left == 0)
		return distSum_right/num_right;
	if(num_right == 0)
		return 100000;
 
	return distSum_right/num_right-distSum_left/num_left;
}
 

bool BestIG_new(Histogram* hist, double* best_t, double* max_gain, double* max_gap){
	
	bool temp = false;
 	
 	TypeCount countLeft;
 	TypeCount countRight;

	list<HistData>::iterator it = hist->begin();
	
 
	double distSum_left = 0;
	double distSum_right = 0;
	
	while(it!=hist->end()){
		countRight.insert(it->type);
		distSum_right = distSum_right+it->value;
		it++;
	}

	countLeft = countRight.defaultCount();


 
	//cout<<"a "<<a2<<" b "<<b2<<" dist "<<distSum_right<<endl;
	
	it = hist->begin();
 
	double totalEntropy = Entropy(&countRight);

	double entropy;
	double gain;
	double dist_before = -1;
	double dist = -1;
 
	while(it!=hist->end()){
 		


		dist_before = dist;
		dist = it->value;
 
		entropy = Entropy(&countLeft, &countRight);
		gain = totalEntropy-entropy;
		
		//cout<<a1<<" "<<b1<<" "<<a2<<" "<<b2<<endl;
		//cout<<dist_before<<" "<<dist<<" "<<gain<<" "<<Gap(distSum_left, distSum_right, a1+b1, a2+b2)<<endl;
		//cout<<endl;
		// if gain increases
		if(gain>*max_gain){
 
			temp = true;
			*max_gain = gain;
			*max_gap = Gap(distSum_left, distSum_right, countLeft.totalCount, countRight.totalCount);
			*best_t = (dist_before+dist)/2;
		}
 
		// if gain ties
		else if(abs(gain-*max_gain)<0.0001){
 
			double gap = Gap(distSum_left, distSum_right, countLeft.totalCount, countRight.totalCount);
			if(gap>*max_gap){

				temp = true;
				*max_gap = gap;
				*best_t = (dist_before+dist)/2;
			}
		}
 
		// change values
		countLeft.add(it->type);
		countRight.subtract(it->type);
 
 
		distSum_right = distSum_right - dist;
		distSum_left = distSum_left + dist;
 
		it++;
 
	}
 
	return temp;
}

bool BestIG(Histogram* hist, double* best_t, double* max_gain, double* max_gap){
	bool temp = false;
 
	list<HistData>::iterator it = hist->begin();
	
	int a1 = 0;
	int b1 = 0;
	int a2 = 0;
	int b2 = 0;
 
	double distSum_left = 0;
	double distSum_right = 0;
	
	while(it!=hist->end()){
		if(it->type==0)
			a2++;
		else
			b2++;
		distSum_right = distSum_right+it->value;
		it++;
	}
 
	//cout<<"a "<<a2<<" b "<<b2<<" dist "<<distSum_right<<endl;
	
	it = hist->begin();
 
	double totalEntropy = Entropy(a2, b2);
	double entropy;
	double gain;
	double dist_before = -1;
	double dist = -1;
 
	while(it!=hist->end()){
 
		dist_before = dist;
		dist = it->value;
 
		entropy = Entropy(a1, b1, a2, b2);
		gain = totalEntropy-entropy;
		//cout<<a1<<" "<<b1<<" "<<a2<<" "<<b2<<endl;
		//cout<<dist_before<<" "<<dist<<" "<<gain<<" "<<Gap(distSum_left, distSum_right, a1+b1, a2+b2)<<endl;
		//cout<<endl;
		// if gain increases
		if(gain>*max_gain){
 
			temp = true;
			*max_gain = gain;
			*max_gap = Gap(distSum_left, distSum_right, a1+b1, a2+b2);
			*best_t = (dist_before+dist)/2;
		}
 
		// if gain ties
		else if(abs(gain-*max_gain)<0.0001){
 
			double gap = Gap(distSum_left, distSum_right, a1+b1, a2+b2);
			if(gap>*max_gap){

				temp = true;
				*max_gap = gap;
				*best_t = (dist_before+dist)/2;
			}
		}
 
		// change values
		if(it->type==0){
			a1++;
			a2--;
		}
		else{
			b1++;
			b2--;
		}
 
 
		distSum_right = distSum_right - dist;
		distSum_left = distSum_left + dist;
 
		it++;
	}
 
	return temp;
}
 

double upperIG(Histogram* hist, double R){
	
	TypeCount countLeft;
 	TypeCount countRight;
 	TypeCount countLeftTemp;
 	TypeCount countRightTemp;

 	//ititialization
	list<HistData>::iterator it = hist->begin();

	while(it!=hist->end()){
		countRight.insert(it->type);
		it++;
	}

	countLeft = countRight.defaultCount();

	it = hist->begin();

	double dist_before = -1;
	double dist = -1;
	double t;
	double entropy;
	double gain;
	double max_gain = 0;
	double totalEntropy = Entropy(&countRight);

	list<HistData>::iterator it_left = it;
	list<HistData>::iterator it_right = it;


	//iteration
	while(it!=hist->end()){
		
		dist_before = dist;
		dist = it->value;
 		t = (dist+dist_before+0.0)/2.0;

		//cout<<t<<endl;

 		countLeftTemp = countLeft;
 		countRightTemp = countRight; 

 		//check left
 		if(it_left!=hist->begin())
 			it_left--;

		while(it_left->value>t-R){	
			
			//cout<<"*  "<<it_left->type<<" "<<(countLeft.count(it_left->type)+0.0)/countLeft.totalCount<<" "<<(countRight.count(it_left->type)+0.0)/countRight.totalCount<<endl;
			if((countLeft.count(it_left->type)+0.0)/countLeft.totalCount<(countRight.count(it_left->type)+0.0)/countRight.totalCount){
				countLeftTemp.subtract(it_left->type);
				countRightTemp.add(it_left->type);
			}

			if(it_left==hist->begin())
				break;

			it_left--;

		}

		//check right
		while(it_right->value<t+R){
			if(it_right==hist->end())
				break;

			//cout<<"** "<<it_right->type<<" "<<(countLeft.count(it_right->type)+0.0)/countLeft.totalCount<<" "<<(countRight.count(it_right->type)+0.0)/countRight.totalCount<<endl;
			if((countLeft.count(it_right->type)+0.0)/countLeft.totalCount>(countRight.count(it_right->type)+0.0)/countRight.totalCount){
				countLeftTemp.add(it_right->type);
				countRightTemp.subtract(it_right->type);
			}

			it_right++;
		}

		entropy = Entropy(&countLeftTemp, &countRightTemp);
		gain = totalEntropy-entropy;

		// //for debugging
		// cout<<" "<<countLeft.count(0)<<" "<<countLeft.count(1)<<endl;
		// cout<<" "<<countRight.count(0)<<" "<<countRight.count(1)<<endl;
		// cout<<" "<<countLeftTemp.count(0)<<" "<<countLeftTemp.count(1)<<endl;
		// cout<<" "<<countRightTemp.count(0)<<" "<<countRightTemp.count(1)<<endl;
		// cout<<gain<<endl;
		// cout<<endl;

		if(gain>max_gain)
			max_gain = gain;

		//move iterators
		countLeft.add(it->type);
		countRight.subtract(it->type);

		it++;
		it_left = it;
		it_right = it;
	}

	return max_gain;
}

bool isDataPure(vector<Data>* dataVector){
	if(dataVector->size()==0)
		return true;
	int type = (*dataVector)[0].type;
	for(int i=0; i<dataVector->size(); i++)
		if((*dataVector)[i].type!=type)
			return false;
	return true;
}

void Shapelet_Discovery(vector<Data>* dataVector, Node* node){
 	
	Stats stats(dataVector);
 
	ShapeletInfo Shapelet;
 	
	vector<double> candidate;

	double max_gain = 0;
	double max_gap = 0;
 	
 	int best_k;
 	int best_l;
 	int best_i;

	int len;

	//for each data
	for(int k=0; k<dataVector->size(); k++){
		Data d = (*dataVector)[k];
		len = d.size();

		for(int l=MINLEN; l<=MAXLEN; l++){
			for(int i=0; i<=len-l; i++){
				candidate.resize(l);
				memcpy(&(candidate[0]), &(d.data[i]), sizeof(double)*l);

				//make Histogram
				Histogram hist;
				for(int t=0; t<dataVector->size(); t++){
					//double distance = sdist(&candidate, &(*dataVector)[t].data);
					double distance = sdist_new(dataVector, &stats, k, t, i, l);
					hist.insert((*dataVector)[t].type, distance);
				}
				if(BestIG_new(&hist, &(Shapelet.best_t), &max_gain, &max_gap)){
					best_k = k;
					best_i = i;
					best_l = l;
					// hist.print();
					// cout<<i<<" "<<j<<" :  gain = "<<max_gain<<", gap = "<<max_gap<<endl;
				}
			}
		}
		cout<<"Round "<<(k<10?"0":"")<<k<<" : t = "<<(Shapelet.best_t)<<", gain = "<<max_gain<<", gap = "<<max_gap<<endl;
	}

	for(int i=0; i<best_l; i++){
		Shapelet.best_S.push_back(((*dataVector)[best_k]).data[i+best_i]);
	}
	node->shapeletInfo = Shapelet;
	// cout<<best_k<<" "<<best_i<<" "<<best_l<<endl;
	// printVector(node->shapeletInfo.best_S);

	//FOR Recursion
	vector <Data> dataVector_left;
	vector <Data> dataVector_right;
	for(int k=0; k<dataVector->size(); k++){

		double distance = sdist(&(node->shapeletInfo.best_S), &((*dataVector)[k].data));
		if(distance<node->shapeletInfo.best_t){
			dataVector_left.push_back((*dataVector)[k]);
		}
		else{
			dataVector_right.push_back((*dataVector)[k]);
		}
	}

	if(isDataPure(&dataVector_left)==false){
		node->left = new Node();
		Shapelet_Discovery(&dataVector_left, node->left);
	}
	else{
		node->left = new Node();
		node->left->type = dataVector_left[0].type;
	}

	if(isDataPure(&dataVector_right)==false){
		node->right = new Node();
		Shapelet_Discovery(&dataVector_right, node->right);
	}
	else{
		node->right = new Node();
		node->right->type = dataVector_right[0].type;
	}
}

void Fast_Shapelet_Discovery(vector<Data>* dataVector, Node* node){
 	
	Stats stats(dataVector);
 
	ShapeletInfo Shapelet;
 
	double max_gain = 0;
	double max_gap = 0;
 	int best_k;
 	int best_l;
 	int best_i;

	vector<double> candidate;
	PruningCandidate pCandidate;

	int len;
	int R;

	//for each data
	for(int k=0; k<dataVector->size(); k++){
		Data d = (*dataVector)[k];
		len = d.size();

		//for every length
		for(int l=MINLEN; l<=MAXLEN; l++){

			pCandidate.clear();

			//for every starting point
			for(int i=0; i<=len-l; i++){
				bool skip = false;

				candidate.resize(l);
				memcpy(&(candidate[0]), &(d.data[i]), sizeof(double)*l);

				//PRUNING!!
				list<PruningCandidateItem>::iterator it = pCandidate.begin();
				while(it!=pCandidate.end()){
					R = sdist_new(dataVector, &stats, k, it->k, i, it->u, l);
					if(upperIG(&(it->hist), R)<0.5*max_gain){
						skip = true;
						break;
					}
					it++;
				}
				if(skip)
					continue;


				//make Histogram

				Histogram hist;
				for(int t=0; t<dataVector->size(); t++){
					//double distance = sdist(&candidate, &(*dataVector)[t].data);
					double distance = sdist_new(dataVector, &stats, k, t, i, l);
					hist.insert((*dataVector)[t].type, distance);
				}
				if(BestIG_new(&hist, &(Shapelet.best_t), &max_gain, &max_gap)){
					best_k = k;
					best_i = i;
					best_l = l;
				}
				pCandidate.insert(&hist, k, i, l);
			}
		}
		cout<<"Round "<<(k<10?"0":"")<<k<<" : t = "<<(Shapelet.best_t)<<", gain = "<<max_gain<<", gap = "<<max_gap<<endl;

	}

	for(int i=0; i<best_l; i++){
		Shapelet.best_S.push_back(((*dataVector)[best_k]).data[i+best_i]);
	}
	node->shapeletInfo = Shapelet;
	// cout<<best_k<<" "<<best_i<<" "<<best_l<<endl;
	// printVector(node->shapeletInfo.best_S);

	//FOR Recursion
	vector <Data> dataVector_left;
	vector <Data> dataVector_right;
	for(int k=0; k<dataVector->size(); k++){

		double distance = sdist(&(node->shapeletInfo.best_S), &((*dataVector)[k].data));
		if(distance<node->shapeletInfo.best_t){
			dataVector_left.push_back((*dataVector)[k]);
		}
		else{
			dataVector_right.push_back((*dataVector)[k]);
		}
	}

	if(isDataPure(&dataVector_left)==false){
		node->left = new Node();
		Fast_Shapelet_Discovery(&dataVector_left, node->left);
	}
	else{
		node->left = new Node();
		node->left->type = dataVector_left[0].type;
	}

	if(isDataPure(&dataVector_right)==false){
		node->right = new Node();
		Fast_Shapelet_Discovery(&dataVector_right, node->right);
	}
	else{
		node->right = new Node();
		node->right->type = dataVector_right[0].type;
	}
}
 
 
/////////////////////////////////////////////////////////////////////////////////////////////////

int setType(Data* data, Node* node){
	if(node->type!=-1){
		return node->type;
	}
	double distance = sdist(&(node->shapeletInfo.best_S), &(data->data));

	if(distance<node->shapeletInfo.best_t){
		return setType(data, node->left);
	}
	else{
		return setType(data, node->right);
	}
		
}

void setTypes(vector<Data>* dataVector, Node* node){
	for(int i=0; i<dataVector->size(); i++){
		Data d = (*dataVector)[i];
		(*dataVector)[i].type2 = setType(&d, node);
		//cout<<setType(&d, node)<<endl;
	}
}

double GetAccuracy(vector <Data>* dataVector2){
	int total = 0;
	int correct = 0;
	for(int k=0; k<dataVector2->size(); k++){
		total++;
		//cout<<k<<" : "<<(*dataVector2)[k].type<<" "<<(*dataVector2)[k].type2<<endl;
		if((*dataVector2)[k].type==(*dataVector2)[k].type2)
			correct++;
	}
	return correct*100.0/total;
}

/////////////////////////////////////////////////////////////////////////////////////////////////


int main(){
	clock_t before;
 
	before = clock();
 	
	Node root;

	vector <Data> dataVector;
 	ifstream inFile;
	inFile.open("gun_train_2", ios::in);
	string s;
 
	while(getline(inFile, s)){
		stringstream ss(s);
		double test;
		int t;
		Data temp;
 
		ss>>t;
		temp.type = t;
		
		while(ss>>test){
			temp.data.push_back(test);
		}
		dataVector.push_back(temp);
	}
 	

	// //FOR DEBUGGING - DATA SET SIZE = 10
	// for(int i=0; i<10; i++){
	// 	getline(inFile, s);
	// 	stringstream ss(s);
	// 	double test;
	// 	int t;
	// 	Data temp;
 
	// 	ss>>t;
	// 	temp.type = t;
		
	// 	while(ss>>test){
	// 		temp.data.push_back(test);
	// 	}
	// 	dataVector.push_back(temp);
	// }
 	
	// for(int i=0; i<dataVector.size(); i++){
	// 	cout<<i<<endl;;
	// 	printVector(dataVector[i].data);
	// 	cout<<endl;
	// }
 	

	// TypeCount typeCount(&dataVector);
 	
 // 	typeCount.print();

	// TypeCount typeCount2 = typeCount.defaultCount();
 	
 // 	typeCount2.print();

 // 	typeCount2.add(1);
 // 	typeCount2.add(2);
 // 	typeCount2.add(2);


 // 	typeCount2.print();

	Fast_Shapelet_Discovery(&dataVector, &root);
 	
 	root.printNodes();

 	//GET TIME
	double result = (double)(clock()-before)/CLOCKS_PER_SEC;
	cout<<"Time = "<<result<<endl;

	/////////TESTING//////////
	
	vector <Data> dataVector2;	

	ifstream inFile2;
	inFile2.open("gun_test", ios::in);
	string s2;

	while(getline(inFile2, s2)){
		stringstream ss(s2);
		double test;
		int t;
		Data temp;

		ss>>t;
		temp.type = t;

		while(ss>>test){
			temp.data.push_back(test);
		}
		dataVector2.push_back(temp);
	}

	setTypes(&dataVector2, &root);

	cout<<"The Accuracy of this Algorithm is : "<<GetAccuracy(&dataVector2)<<"%"<<endl;


	//////////////////////upperIG TESTING////////////////////////////

	// Histogram hist;
	// hist.insert(0, 0);
	// hist.insert(1, 1);
	// hist.insert(0, 2);
	// hist.insert(0, 3);
	// hist.insert(0, 4);
	// hist.insert(0, 5);
	// hist.insert(0, 6);
	// hist.insert(0, 7);

	// cout<<upperIG(&hist, 1);

	/////////////////////////////////////////////////////////////////

	//////////////////////sdist_new and dist TESTING///////////////////////////
 
	// Data data1;
	// Data data2;
 
	// data1.data.push_back(1);
	// data1.data.push_back(2);
	// data1.data.push_back(3);
 
	// data1.data.push_back(4);
	// data1.data.push_back(5);
	// data1.data.push_back(4);
	// data1.data.push_back(5);
 
	// data1.data.push_back(6);
	// data1.data.push_back(7);
	// data1.data.push_back(8);
 
	// data2.data.push_back(9);
	// data2.data.push_back(8);
	// data2.data.push_back(6);
	// data2.data.push_back(7);
	// data2.data.push_back(4);
	// data2.data.push_back(3);
	// data2.data.push_back(5);
	// data2.data.push_back(2);
	// data2.data.push_back(0);
	// data2.data.push_back(1);
 
 
	// Data data3;
	// data3.data.push_back(4);
	// data3.data.push_back(5);
	// data3.data.push_back(4);
	// data3.data.push_back(5);
 
	// vector<Data> dataVector2;
	// dataVector2.push_back(data1);
	// dataVector2.push_back(data2);
 
	// Stats stats2(&dataVector2);
 
	// cout<<sdist(&data3.data, &data2.data)<<endl;
	// cout<<endl;
	// cout<<sdist_new(&dataVector2, &stats2, 0, 1, 3, 4);
 
	/////////////////////////////////////////////////////////////////////////////
 
 
 
	// Histogram hist_temp;
	// hist_temp.insert(1, 1.321);
	// hist_temp.insert(1, 2.321);
	// hist_temp.insert(1, 2.564);
	// hist_temp.insert(1, 3.124);
	// hist_temp.insert(0, 5.123);
	// hist_temp.insert(1, 6.123);
	// hist_temp.insert(0, 6.431);
	// hist_temp.insert(0, 7.123);
	// hist_temp.insert(0, 8.231);
	// hist_temp.insert(0, 9.231);
 
	// hist_temp.print();
	// double best_t = 0;
	// double max_gain = 0;
	// double max_gap = 0;
	// BestIG(&hist_temp, &best_t, &max_gain, &max_gap);
 
	// cout<<endl;
 
	// cout<<best_t<<" "<<max_gain<<" "<<max_gap<<endl;
	// vector<double> test;
	// test.push_back(1);
	// test.push_back(2);
	// test.push_back(3);
	// test.push_back(2);
	// test.push_back(4);
	// test.push_back(1);
	// test.push_back(2);
	// test.push_back(3);
	// test.push_back(1);
	// test.push_back(4);
 
 
 
	// ShapeletInfo info;
	// info.best_S = test;
	// info.best_t = 10;
	// printShapelet(info);
 
 
 
	// // vector<double> test2 = zNorm(&test, 0, test.size());
	// // printVector(test2);
	
 
	// vector<double> test2;
	// test2.push_back(3);
	// test2.push_back(2);
	// test2.push_back(4);
	// test2.push_back(1);
 
	// cout<<sdist(&test2, &test);
}