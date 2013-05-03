#include "stdafx.h"
#include "vector3d.h"

vector3d::vector3d(){
	x = 0;
	y = 0;
	z = 0;
}

vector3d::vector3d(double vx, double vy, double vz){
	x = vx;
	y = vy;
	z = vz;
}

vector3d& vector3d::operator = (vector3d& v){
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
}

double vector3d::getNorm(){
	return sqrt(x*x + y*y + z*z);
}