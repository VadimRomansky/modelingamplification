#include "math.h"

#include "vector3d.h"

Vector3d::Vector3d(){
	x = 0;
	y = 0;
	z = 0;
}

Vector3d::Vector3d(double vx, double vy, double vz){
	x = vx;
	y = vy;
	z = vz;
}

double Vector3d::getNorm(){
	return sqrt(x*x + y*y + z*z);
}

Vector3d Vector3d::operator+(const Vector3d& vector){
	return Vector3d(x + vector.x, y + vector.y, z + vector.z);
}

Vector3d Vector3d::operator-(const Vector3d& vector){
	return Vector3d(x - vector.x, y - vector.y, z - vector.z);
}

Vector3d Vector3d::operator*(const double& value){
	return Vector3d(x*value, y*value, z*value);
}

Vector3d Vector3d::operator/(const double& value){
	return Vector3d(x/value, y/value, z/value);
}

double Vector3d::scalarMult(const Vector3d& vector){
	return x*vector.x + y*vector.y + z*vector.z;
}

Vector3d Vector3d::vectorMult(const Vector3d& vector){
	double x = y*vector.z - z*vector.y;
	double y = -x*vector.z + z*vector.x;
	double z = x*vector.y - y*vector.x;

	return Vector3d(x, y, z);
}