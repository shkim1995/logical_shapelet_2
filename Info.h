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


int MAXLEN = 60;
int MINLEN = 30;
int CANLEN = 1;

 
typedef struct _ShapeletInfo
{
	vector<double> best_S;
	double best_t;
	double maxGain;
	double maxGap;
 
}ShapeletInfo;