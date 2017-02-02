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


class Stats{
public:
 
	double** S; //array of SumArray
	double** S2; // array of SquareSumArray
	double**** M;
 
	Stats(vector<Data>* dataVector){
		
		//Inititate S, S2
		int dataLen = dataVector->size();
		S = new double* [dataLen];
		S2 = new double* [dataLen];
 
		int len;
		double sum;
		double squreSum;
		double temp;
 
		for(int i=0; i<dataLen; i++){
			len = (*dataVector)[i].size();
			sum = 0;
			squreSum = 0;
 
			S[i] = new double[len+1];
			S2[i] = new double[len+1];
			S[i][0] = 0;
			S2[i][0] = 0;
 
			for(int j=0; j<len; j++){
				temp = ((*dataVector)[i].data)[j];
				sum = sum + temp;
				squreSum = squreSum + temp*temp;
				S[i][j+1] = sum;
				S2[i][j+1] = squreSum;
			}
		}
 
		//Inititate M
		int x = 0;
		int y = 0;
 
		M = new double*** [dataLen];
		for(int i=0; i<dataLen; i++)
			M[i] = new double** [dataLen];
 
		for(int i=0; i<dataLen; i++){
			for(int j=0; j<dataLen; j++){
				
				x = (*dataVector)[i].size();
				y = (*dataVector)[j].size();
 
				//inititalization
				M[i][j] = new double* [x+1];
				for(int u=0; u<x+1; u++)
					M[i][j][u] = new double[y+1];
 
				//boundary condition
				for(int u=0; u<x+1; u++)
					M[i][j][u][0] = 0;
				for(int v=0; v<y+1; v++)
					M[i][j][0][v] = 0;
 
				for(int v=1; v<y+1; v++){
					for(int u=v; u<x+1; u++){
						M[i][j][u][v] = M[i][j][u-1][v-1]+((*dataVector)[i].data)[u-1]*((*dataVector)[j].data)[v-1];
						
					}
				}
 
				for(int u=1; u<x+1; u++){
					for(int v=u; v<y+1; v++){
						M[i][j][u][v] = M[i][j][u-1][v-1]+((*dataVector)[i].data)[u-1]*((*dataVector)[j].data)[v-1];
						
					}
				}
			}
		}
 
	}
 
};