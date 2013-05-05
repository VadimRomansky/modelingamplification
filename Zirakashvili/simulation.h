#ifndef SIMULATION_H
#define SIMULATION_H

#include "stdafx.h"

class Simulation{
public:
	double timeStep;
	int iterationNumber;
	int particlesNumber;
	int allParticlesNumber;
	double U0;
	double Rtot;
	double Rloc;
	double density0;
	double B0;
	double upstreamR;
	double downstreamR;
	double deltaR;
	double deltaTheta;
	double deltaPhi;
	double temperature;
	double gridParameter;
	double epsilonR;
	double simulationTime;
	double deltaT;
	bool kolmogorovCascading;
	bool resonantInstability;
	bool bellInstability;
	int A;
	int Z;

	bool smallAngleScattering;

	double energy;
	double theorEnergy;

	double momentumZ;
	double theorMomentumZ;
	double momentumX;
	double theorMomentumX;
	double momentumY;
	double theorMomentumY;
	double particlesWeight;

	int rgridNumber;
	int shockWavePoint;
	int currentShockWavePoint;

	Particle* Particles;
	SpaceBin** bins;
	double** pressureSpectralDensity;
	std::list <Particle*> startPDF;
	std::vector<Particle*> introducedParticles;

	double* averageVelocity;

	Simulation();
	~Simulation();

	void initializeProfile();
	void simulate();
};

#endif