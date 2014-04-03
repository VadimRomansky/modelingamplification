#include "math.h"
#include "constants.h"

double u(double x){
	double u1, u2;

	//u1 = 1;
	u1 = 3.0E5;
	u2 = 0.25*u1;

	if  (x < 0) {
		return u1;
	} else {
		return u2;
	}
}

double kappa(double x,double y){
	//double kappa0 = 0.1;
	double B0 = 6E-6;
	double kappa0 = speed_of_light*speed_of_light/(electron_charge*B0);
    return kappa0*exp(y);
}

double Qinj(double x, double dx, double y,double dy){
	if (x < 0.1*dx && x > -0.1*dx && y < 0.5*dy && y > -0.5*dy) {
		double middleDensity = 0.0000000000000000000000016;
		return 0.1*middleDensity*abs(u(x))/(exp(3*y)*massProton*dy);;
	} else {
		return 0;
	}
}

double lim(double fp, double f0, double fm){
	double r;

	if  ((fp-f0) == 0.0) {
		return 2;
	} else {
		r= (f0 - fm)/(fp - f0);
		return (r + abs(r))/(1 + abs(r));
	}
}