#include "math.h"

#include "vector3d.h"

Vector3d::Vector3d(double vx, double vy, double vz){
	x = vx;
	y = vy;
	z = vz;
}

double Vector3d::getNorm(){
	return sqrt(x*x + y*y + z*z);
}

Vector3d Vector3d::operator +(const Vector3d& vector){
	return Vector3d(x + vector.x, y + vector.y, z + vector.z);
}

Vector3d Vector3d::operator -(const Vector3d& vector){
	return Vector3d(x - vector.x, y - vector.y, z - vector.z);
}

double Vector3d::scalarMult(const Vector3d& vector){
	return x*vector.x + y*vector.y + z*vector.z;
}