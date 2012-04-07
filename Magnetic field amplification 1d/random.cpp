#include "stdafx.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "random.h"

double randomMaxwell(double temperature, double mass){
	double summ = 0;
	for(int i = 0; i < randomSeed; i++){
	    summ = summ + (rand()%randomSeed - randomSeed/2)*12*sqrt(mass*kBoltzman*temperature)/(randomSeed);
	}
	summ = summ/sqrt(randomSeed*1.0);
    return summ;
	//return randomGauss(0.0, sqrt(2*mass*kBoltzman*temperature));
}

double randomGauss(double a, double sigma){
	double x = uniRandom();
	double y = uniRandom();
	double t = x*x + y*y;
	while((t > 1) || (abs(t) < 100000000000000*DBL_EPSILON)){
		x = uniRandom();
		t = x*x + y*y;
		if(( t <= 1) && (t > 100000000000000*DBL_EPSILON)){
			break;
		}
		y = uniRandom();
		t = x*x + y*y;
	}
	return a + sigma*sqrt(-2*log(t)/t);
}

double uniRandom(){
	//srand ( time(NULL) );
	double a = (rand()%randomSeed);
	//double a = randomSeed/2;
	a = a/randomSeed;
	return a;
	//return (rand()%randomSeed)/randomSeed;
}