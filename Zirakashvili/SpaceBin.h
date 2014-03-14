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

	double xi;

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
	int numberPhi;
	int numberTheta;

	bool smallAngleScattering;

	std::list <Particle*> detectedParticlesR2;
	std::list <Particle*> detectedParticlesR1;
	std::list <Particle*> detectedParticlesTheta2;
	std::list <Particle*> detectedParticlesTheta1;
	std::list <Particle*> detectedParticlesPhi2;
	std::list <Particle*> detectedParticlesPhi1;

	std::list <Particle*> particles;
	double* magneticField;
	std::list <Particle*>* sortedParticles;

	double minK;
	double maxK;
	double averageVelocity;

	double crFlux;

	double centralMomentum;

	SpaceBin();
	SpaceBin(double r, double xi, double u, double rho, double p);
	~SpaceBin();
	int propagateParticle(Particle* particle ,double& time, double timeStep, const int rgridNumber);
	double getFreePath(Particle* particle);
	void makeOneStep(Particle* particle, double deltat, double& time);
	static int binByCoordinates(double r, double theta, double phi, double r0, double deltar, double deltatheta, double deltaphi, const int rgridNumber); 
	void scattering(Particle* particle, double maxTheta);
	void multiplyParticleWeights(double value);
	void resetDetectors();

	void resetProfile(double massFlux0, double momentaFlux0, double energyFlux0,double density0, double U0);
	double pressureSpectralDensity(int j, double density0, double U0, int Z, int A);
	double getMaxP();
	double getMinP();
	double getMaxLocalP();
	double getMinLocalP();

	double getEnergy();

	void updateTemperature(double* distribution, double deltap);

	void largeAngleScattering(Particle* particle, double& time, double timeStep);

	//void operator=(Xbin bin);

};


#endif