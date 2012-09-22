#include "stdafx.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "random.h"

double randomMaxwell(double temperature, double mass){
	//srand ( time(NULL) );
	////TODO найти ошибку!
	/*double summ = 0;
	for(int i = 0; i < randomSeed; i++){
	    summ = summ + (rand()%randomSeed - randomSeed/2)*12*sqrt(mass*kBoltzman*temperature)/(randomSeed);
	}
	summ = summ/sqrt(randomSeed*1.0);

	return summ;*/

	double x = (uniRandom() - 0.5)*2;
	double y = (uniRandom() - 0.5)*2;
	double s = sqrt(x*x+y*y);
	while((s == 0) || (s > 1)){
		x = (uniRandom() - 0.5)*2;
		y = (uniRandom() - 0.5)*2;
		s = sqrt(x*x+y*y);
	}

	double xi = x*sqrt(-2*log(s)/s);


    return xi*sqrt(mass*kBoltzman*temperature);
}

double uniRandom(){
	//srand ( time(NULL) );
	double a = (rand()%randomSeed);
	a = a/randomSeed;
	return a + 1.0/(2.0*randomSeed);
}