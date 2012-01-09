#include "stdafx.h"
#include "vector3d.h"

vector3d::vector3d(double vx, double vy, double vz){
	x = vx;
	y = vy;
	z = vz;
}

double vector3d::getNorm() const{
	return sqrt(x*x + y*y + z*z);
}