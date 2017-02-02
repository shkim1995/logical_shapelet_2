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