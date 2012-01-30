#ifndef XBIN_H
#define XBIN_H

#include "stdafx.h"
#include "constants.h"
#include "random.h"
#include "BinFlux.h"
//#include "particle2.h"
#include <list>

class Particle;

class SpaceBin{
public:
	static const int defaultTimeRelation = 5;

	double r;
	double r1;
	double r2;
	double theta;
	double theta1;
	double theta2;
	double phi;
	double phi1;
	double phi2;

	double U;
	double UTheta;
	double UPhi;
	double density;
	double volume;
	double temperature;
	double B;
	double B0;
	double W;

	double initialMomentum;

	int numberR;
	int numberPhi;
	int numberTheta;

	BinFlux crMassFlux;
	BinFlux particleMassFlux;
	BinFlux bulkMassFlux;
	BinFlux particleMomentaFlux;
	BinFlux bulkMomentaFlux;
	BinFlux thMomentaFlux;
	BinFlux magMomentaFlux;
	BinFlux particleEnergyFlux;
	BinFlux bulkEnergyFlux;
	BinFlux thEnergyFlux;
	BinFlux magEnergyFlux;

	BinFlux massFlux;
	BinFlux momentaFlux;
	BinFlux energyFlux;

	std::list <Particle*> detectedParticlesR2;
	std::list <Particle*> detectedParticlesR1;
	std::list <Particle*> detectedParticlesTheta2;
	std::list <Particle*> detectedParticlesTheta1;
	std::list <Particle*> detectedParticlesPhi2;
	std::list <Particle*> detectedParticlesPhi1;
	double* magneticField;
	std::list <Particle*>* sortedParticles;

	double minK;
	double maxK;

	SpaceBin();
	SpaceBin(double r, double theta, double phi, double deltar, double deltatheta, double deltaphi, double u, double rho, double utheta,double uphi, double temperature, double b, int i, int j, int k);
	~SpaceBin();
	int* propagateParticle(Particle* particle ,double& time, double timeStep);
	double getFreePath(Particle* particle);
	void makeOneStep(Particle* particle, double deltat, double& time);
	bool isInBin(Particle* particle);
	bool isInThisOrNear(double r, double theta, double phi);
	static int* binByCoordinates(double r, double theta, double phi, double r0, double deltar, double deltatheta, double deltaphi); 
	void scattering(Particle* particle, double maxTheta);
	void updateFluxes();
	void updateCosmicRayFluxes();
	void updateBulkFluxes();
	void updateThermalFluxes();
	void updateMagneticFluxes();
	void updateMagneticField();
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
	void detectParticleTheta1(Particle* particle);
	void detectParticleTheta2(Particle* particle);
	void detectParticlePhi1(Particle* particle);
	void detectParticlePhi2(Particle* particle);

	//void operator=(Xbin bin);

};


#endif