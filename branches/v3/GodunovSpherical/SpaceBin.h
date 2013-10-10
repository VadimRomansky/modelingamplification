#ifndef XBIN_H
#define XBIN_H

#include "stdafx.h"
#include "constants.h"
#include "random.h"
#include "matrix3d.h"
#include <list>

class Particle;

class SpaceBin{
public:
	static const int defaultTimeRelation = 5;

	double r;
	double r1;
	double r2;

	double U;
	double density;
	double volume;
	double temperature;
	double pressure;

	double B0;

	double u1;
	double u2;
	double u3;

	double initialMomentum;

	Matrix3d* matrix;
	Matrix3d* invertMatrix;

	int numberR;

	bool smallAngleScattering;

	SpaceBin();
	SpaceBin(double r, double deltaR, double u, double rho, double temperature, double B0, int i, bool scattering);
	~SpaceBin();
	int propagateParticle(Particle* particle ,double& time, double timeStep, const int rgridNumber);
	double getFreePath(Particle* particle);
	void makeOneStep(Particle* particle, double deltat, double& time);
	static int binByCoordinates(double r, double theta, double phi, double r0, double deltar, double deltatheta, double deltaphi, const int rgridNumber); 
	void scattering(Particle* particle, double maxTheta);
	void multiplyParticleWeights(double value);
	void resetDetectors();

	double getEnergy();

	void largeAngleScattering(Particle* particle, double& time, double timeStep);

};


#endif