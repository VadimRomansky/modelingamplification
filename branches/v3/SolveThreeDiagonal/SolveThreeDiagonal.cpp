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
	double alpha = f[0]/middle[0];
	double betta = - upper[0]/middle[0];
	x[0] = 0;
	for(int i = 1; i < rgridNumber - 1; ++i){
		x[i] = 0;
		alpha = (f[i] - lower[i-1]*alpha)/(lower[i-1]*betta + middle[i]);
		alertNaNOrInfinity(alpha, "aplpha = NaN");
		betta = - upper[i]/(lower[i-1]*betta + middle[i]);
		alertNaNOrInfinity(betta, "betta = NaN");
	}
	x[rgridNumber - 1] = (f[rgridNumber - 1] - lower[rgridNumber - 2]*alpha)/(lower[rgridNumber-2]*betta + middle[rgridNumber - 1]);
	alertNaNOrInfinity(x[rgridNumber - 1], "x[rgridNumber - 1] = NaN");
	x[rgridNumber - 2] = (f[rgridNumber - 1] - middle[rgridNumber - 1]*x[rgridNumber - 1])/(lower[rgridNumber - 2]);
	alertNaNOrInfinity(x[rgridNumber - 2], "x[rgridNumber - 2] = NaN");
	for(int i = rgridNumber - 3; i >= 0; --i){
		x[i] = (f[i + 1] - middle[i + 1]*x[i + 1] - upper[i + 1]*x[i+2])/lower[i];
		alertNaNOrInfinity(x[i],"x = NaN");
		/*if(abs(x[i]) > epsilon){
			printf("aaaaaa\n");
		}*/
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	double* upper = new double[rgridNumber - 1];
	double* middle = new double[rgridNumber];
	double* lower = new double[rgridNumber];
	double* f = new double[rgridNumber];
	double* x = new double[rgridNumber];
	
	upper[0] = 2;
	upper[1] = 3;
	upper[2] = 1;
	middle[0] = 7;
	middle[1] = 1;
	middle[2] = 1;
	middle[3] = 1;
	lower[0] = 1;
	lower[1] = 2;
	lower[2] = 1;

	f[0] = 1;
	f[1] = 0;
	f[2] = 0;
	f[3] = 0;

	solveThreeDiagonal(middle, upper, lower, f, x);

	for(int i = 0; i < rgridNumber; ++i){
		printf("%lf\n", x[i]);
	}

	delete[] upper;
	delete[] middle;
	delete[] lower;
	delete[] f;
	delete[] x;
	return 0;
}

