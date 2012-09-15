#include "stdafx.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "random.h"

double randomMaxwell(double temperature, double mass){
	//srand ( time(NULL) );
	////TODO найти ошибку!
	double summ = 0;
	for(int i = 0; i < randomSeed; i++){
	    summ = summ + (rand()%randomSeed - randomSeed/2)*12*sqrt(mass*kBoltzman*temperature)/(randomSeed);
	}
	summ = summ/sqrt(randomSeed*1.0);
    return summ;
}

double uniRandom(){
	//srand ( time(NULL) );
	double a = (rand()%randomSeed);
	a = a/randomSeed;
	//return a;
	return a + 1.0/(2.0*randomSeed);
}