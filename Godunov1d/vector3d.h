#ifndef VECTOR3D_H
#define VECTOR3D_H

class vector3d{
public:
	double x;
	double y;
	double z;
	vector3d();
	vector3d(double vx, double vy, double vz);
	vector3d& operator = (vector3d& v);
	double getNorm();
};

#endif