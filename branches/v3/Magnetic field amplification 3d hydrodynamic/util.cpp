#include "stdafx.h"
#include "util.h"
#include "vector3d.h"
#include "constants.h"

double power(double v, double p){
	return exp(p*log(v));
}

double sqr(double v){
	return v*v;
}

int lowerInt(double v){
	int i = (int) v;

	if(1.0*i >= v){
		i = i - 1;
	}
	return i;
}

vector3d summVelocity(vector3d v, double u){
	double c2 = speed_of_light*speed_of_light;
	double vz = (v.z - u)/(1 - v.z*u/c2);
	double vy = v.y*sqrt(1 - u*u/c2)/(1 - v.z*u/c2);
	double vx = v.x*sqrt(1 - u*u/c2)/(1 - v.z*u/c2);

	return vector3d(vx,vy,vz);
}

bool order(double a, double b, double c){
	return ((a <= b) && (b <= c));
}

double angleDelta(double phi1, double phi2){
	double delta1 = phi2 - phi1;
	double delta2 = delta1 + 2*pi;
	double delta3 = delta1 - 2*pi;
	if(abs(delta1) < abs(delta2)){
		if(abs(delta1) < abs(delta3)){
			return delta1;
		} else {
			return delta3;
		}
	} else {
		if(abs(delta2) < abs(delta3)){
			return delta2;
		} else {
			return delta3;
		}
	}
}

double max(double a, double b){
	if(a >= b){
		return a;
	} else {
		return b;
	}
}

double min4(double a, double b,double c, double d){
	if( a <  b ){
		if( c < d){
			if( a < c){
				return a;
			} else {
				return c;
			}
		} else {
			if( a < d){
				return a;
			} else {
				return d;
			}
		}
	} else {
		if( c < d){
			if( b < c){
				return b;
			} else {
				return c;
			}
		} else {
			if( b < d){
				return b;
			} else {
				return d;
			}
		}
	}
}

double maxwell(double momentum, double mass, double temperature){
	return exp(-momentum*momentum/(2*mass*kBoltzman*temperature))/sqrt(2*pi*mass*kBoltzman*temperature);
}

void alertNaNOrInfinity(double value, const char* s){
	if(value != value || 0*value != 0*value){
		printf(s);
		printf("\n");
	}
}