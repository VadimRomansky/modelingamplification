#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "stdlib.h"
#include "vector"

#include "particle.h"

class Simulation{
	int xnumber;
	int ynumber;
	int znumber;
	int particlesNumber;

	double deltaX;
	double deltaY;
	double deltaZ;

	double* xgrid;
	double* ygrid;
	double* zgrid;
	double* middleXgrid;
	double* middleYgrid;
	double* middleZgrid;


	Vector3d*** Efield;
	Vector3d*** Bfield;

	Vector3d*** newEfield;
	Vector3d*** newBfield;

	Vector3d*** tempEfield;
	Vector3d*** tempBfield;

	std::vector<Particle*> particles;

	Simulation();
	~Simulation();

	void initialize();
	void createArrays();
	Vector3d correlationTempEfield(Particle* particle);
	Vector3d correlationBfield(Particle* particle);
	
	Vector3d correlationField(Particle* particle, Vector3d*** field);
	Vector3d correlationFieldWithBin(Particle* particle, Vector3d*** field, int i, int j, int k);
	double correlationBspline(const double& x, const double&  dx, const double& leftx, const double& rightx);
};

#endif