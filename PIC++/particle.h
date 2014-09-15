#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "vector3d.h"

enum ParticleTypes {PROTON, ELECTRON};

class Particle{
public:
	double mass;
	double charge;
	double weight;
	ParticleTypes type;

	Vector3d coordinates;

	Vector3d momentum;

	Matrix3d rotationTensor;


	double dx;
	double dy;
	double dz;

	Particle(double m, double q, double w, ParticleTypes type, double x0, double y0, double z0, double px0, double py0, double pz0, double dx0, double dy0, double dz0);
	Particle(const Particle& particle);

	double shapeFunctionX(const double& xvalue);
	double shapeFunctionY(const double& yvalue);
	double shapeFunctionZ(const double& zvalue);

	double shapeFunction(const double& xvalue, const double& yvalue, const double& zvalue);

	double momentumAbs();
	Vector3d velocity(double c);
	double velocityX(double c);
	double velocityY(double c);
	double velocityZ(double c);
	double gammaFactor(double c);

	void setMomentumByV(Vector3d v, double c);
};

#endif