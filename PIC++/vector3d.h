#ifndef VECTOR3D_H
#define VECTOR3D_H

class Vector3d{
public:
	double x;
	double y;
	double z;
	Vector3d(double vx, double vy, double vz);
	double getNorm();
	Vector3d operator-(Vector3d vector);
	Vector3d operator+(Vector3d vector);
};

#endif