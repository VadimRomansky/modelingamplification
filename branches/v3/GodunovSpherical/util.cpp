#include "stdafx.h"
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

double max(double a, double b){
	if(a >= b){
		return a;
	} else {
		return b;
	}
}

double min(double a, double b){
	if(a >= b){
		return b;
	} else {
		return a;
	}
}

void alertNaNOrInfinity(double value, const char* s){
	if(value != value || 0*value != 0*value){
		printf(s);
		printf("\n");
	}
}

void alertNotPositive(double value, const char* s){
	if(value <= 0){
		printf(s);
		printf("\n");
	}
}