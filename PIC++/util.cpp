#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "constants.h"

double power(double v, double p){
	return exp(p*log(v));
}

double sqr(double v){
	return v*v;
}

double cube(double v){
	return v*v*v;
}

double max2(double a, double b){
	if(a >= b){
		return a;
	} else {
		return b;
	}
}

double min2(double a, double b){
	if(a >= b){
		return b;
	} else {
		return a;
	}
}

void alertNaNOrInfinity(double value, const char* s){
	if(value != value || 0*value != 0*value){
		printf("%s", s);
		printf("\n");
		exit(0);
		//Sleep(1000);
	}
}

void alertNotPositive(double value, const char* s){
	if(value <= 0){
		printf("%s", s);
		printf("\n");
		exit(0);
	}
}

void alertNegative(double value, const char* s){
	if(value < 0){
		printf("%s", s);
		printf("\n");
		exit(0);
	}
}

void solveSpecialMatrix(double** const leftHalf, double* const rightPart, double* const output){
	for(int j = 0; j < 2; ++j){
		leftHalf[2][j] /= leftHalf[2][2];
	}
	rightPart[2] /= leftHalf[2][2];
	leftHalf[2][2] = 1;

	for(int i = 0; i < 2; ++i){
		for(int j = 0; j < 2; ++j){
			leftHalf[i][j] -= leftHalf[2][j]*leftHalf[i][2];
		}
		rightPart[i] -= rightPart[2]*leftHalf[i][2];
		leftHalf[i][2] = 0;
	}

	leftHalf[1][0] /= leftHalf[1][1];
	rightPart[1] /= leftHalf[1][1];
	leftHalf[1][1] = 1;

	leftHalf[0][0] -= leftHalf[1][0]*leftHalf[0][1];
	rightPart[0] -= rightPart[1]*leftHalf[0][1];
	leftHalf[0][1] = 0;

	rightPart[0] /= leftHalf[0][0];

	leftHalf[0][0] = 1;

	output[0] = rightPart[0];
	for(int i = 1; i < 6; ++i){
		output[i] = rightPart[i];
		for(int j = 0; j < min2(i,3); ++j){
			output[i] -= leftHalf[i][j]*output[j];
		}
	}
}

double coordinateDifference(double* const a, double* const b, double dt){
	double result = 0;
	for(int i = 0; i < 3; ++i){
		result += fabs(a[i] - b[i]);
	}

	for(int i = 3; i < 6; ++i){
		result += fabs(a[i] - b[i])*dt;
	}
	return result;
}
