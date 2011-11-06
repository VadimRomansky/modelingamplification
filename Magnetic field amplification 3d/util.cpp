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

int trunc(double v){
	int i = (int) v;

	if(1.0*i > v){
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