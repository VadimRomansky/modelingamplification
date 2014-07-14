#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "vector3d.h"

class Particle{
public:
	double mass;
	double charge;

	Vector3d coordinates;

	Vector3d momentum;

	double dx;
	double dy;
	double dz;

	Particle(double m, double q, double x0, double y0, double z0, double px0, double py0, double pz0, double dx0, double dy0, double dz0);

	double shapeFunctionX(const double& xvalue);
	double shapeFunctionY(const double& yvalue);
	double shapeFunctionZ(const double& zvalue);

	double shapeFunction(const double& xvalue, const double& yvalue, const double& zvalue);

	double momentumAbs();
	double velocityX();
	double velocityY();
	double velocityZ();
};

#endif