#ifndef XBIN_H
#define XBIN_H

#include "stdafx.h"
#include "constants.h"
#include "random.h"
#include "BinFlux.h"
#include "matrix3d.h"
//#include "particle2.h"
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
	double B;
	double B0;
	double W;

	double initialMomentum;

	int numberR;

	bool smallAngleScattering;


	std::list <Particle*> detectedParticlesR2;
	std::list <Particle*> detectedParticlesR1;

	std::list <Particle*> particles;
	double* magneticField;
	std::list <Particle*>* sortedParticles;

	std::vector<double> particleMomentaZ;
	std::vector<double> particleWeights;

	double minK;
	double maxK;
	double averageVelocity;

	double crFlux;

	double centralMomentum;

	SpaceBin();
	SpaceBin(double r, double deltar, double u, double rho, double temperature, double b, int i, bool scattering);
	~SpaceBin();
	int propagateParticle(Particle* particle ,double& time, double timeStep, const int rgridNumber);
	//double getFreePath(Particle* particle);
	double getFreeTime(Particle* particle);
	void makeOneStep(Particle* particle, double deltat, double& time);
	bool isInBin(Particle* particle);
	bool isInThisOrNear(double r);
	static int binByCoordinates(double r, double r0, double deltar, const int rgridNumber); 
	void scattering(Particle* particle, double maxTheta);

	void multiplyParticleWeights(double value);
	void sortParticles(double minK, double maxK);
	void resetDetectors();

	void resetProfile(double massFlux0, double momentaFlux0, double energyFlux0,double density0, double U0);
	double pressureSpectralDensity(int j, double density0, double U0, int Z, int A);
	double getMaxP();
	double getMinP();
	double getMaxLocalP();
	double getMinLocalP();

	void detectParticleR1(Particle* particle);
	void detectParticleR2(Particle* particle);

	void updateCosmicRayBoundMomentum();
	void updateTemperature(double* distribution, double deltap);

	void largeAngleScattering(Particle* particle, double& time, double timeStep);

	//void operator=(Xbin bin);

};


#endif