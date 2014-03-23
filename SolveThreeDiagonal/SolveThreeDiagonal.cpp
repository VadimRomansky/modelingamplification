// SolveThreeDiagonal.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

const int rgridNumber = 4;

void alertNaNOrInfinity(double value, const char* s){
	if(value != value || 0*value != 0*value){
		printf(s);
		printf("\n");
	}
}

void solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x){
	double* alpha = new double[rgridNumber];
	double* beta = new double[rgridNumber];
	alpha[1] = -upper[0]/middle[0];
	beta[1] = f[0]/middle[0];
	for(int i = 2; i < rgridNumber; ++i){
		double temp = lower[i-1]*alpha[i-1] + middle[i-1];
		alpha[i] = -upper[i-1]/temp;
		beta[i] = (f[i-1] - lower[i-1]*beta[i-1])/temp;
	}

	x[rgridNumber - 1] = (f[rgridNumber-1] - lower[rgridNumber-1]*beta[rgridNumber-1])/(lower[rgridNumber-1]*alpha[rgridNumber-1] + middle[rgridNumber-1]);
	alertNaNOrInfinity(x[rgridNumber-1],"x = NaN");

	for(int i = rgridNumber - 2; i >= 0; --i){
		x[i] = alpha[i+1]*x[i+1] + beta[i+1];
		alertNaNOrInfinity(x[i],"x = NaN");
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	double* upper = new double[rgridNumber];
	double* middle = new double[rgridNumber];
	double* lower = new double[rgridNumber];
	double* f = new double[rgridNumber];
	double* x = new double[rgridNumber];

	/*
	
	7 2 0 0     1
	1 1 3 0     0
	0 2 1 1     0
	0 0 1 1     0
	

	x 
	=

	0.1426
	0
	-0.47619
	0.47619
	*/
	
	upper[0] = 2;
	upper[1] = 3;
	upper[2] = 1;
	middle[0] = 7;
	middle[1] = 1;
	middle[2] = 1;
	middle[3] = 1;
	lower[1] = 1;
	lower[2] = 2;
	lower[3] = 1;

	f[0] = 1;
	f[1] = 0;
	f[2] = 0;
	f[3] = 0;

	solveThreeDiagonal(middle, upper, lower, f, x);

	for(int i = 0; i < rgridNumber; ++i){
		printf("%lf\n", x[i]);
	}
	printf("\n");


	//solving matrix for approximations three points like parabol
	/*
	4 2 1   1
	9 3 1   2
	16 4 1  3

	x = 
	0
	1
	-1
	*/

	double p1 = 2;
	double p2 = 3;
	double p3 = 4;
	double y1 = 1;
	double y2 = 2;
	double y3 = 3;

	double det = p1*p1*p2 + p2*p2*p3 + p3*p3*p1 - p3*p3*p2 - p1*p1*p3 - p2*p2*p1;
	double detA = y1*p2 + y2*p3 + y3*p1 - y3*p2 - y1*p3 - y2*p1;
	double detB = p1*p1*y2 + p2*p2*y3 + p3*p3*y1 - p3*p3*y2 - p1*p1*y3 - p2*p2*y1;
	double detC = p1*p1*p2*y3 + p2*p2*p3*y1 + p3*p3*p1*y2 - p3*p3*p2*y1 - p1*p1*p3*y2 - p2*p2*p1*y3;
	
	double a = detA/det;
	double b = detB/det;
	double c = detC/det;

	printf("a = %lf\nb = %lfc = %lf\n",a,b,c);

	delete[] upper;
	delete[] middle;
	delete[] lower;
	delete[] f;
	delete[] x;
	return 0;
}

