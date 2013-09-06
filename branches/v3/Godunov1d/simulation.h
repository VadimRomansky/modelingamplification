#ifndef SIMULATION_H
#define SIMULATION_H

#include "stdafx.h"
#include "SpaceBin.h"
#include "particle.h"

class Simulation{
public:
	double timeStep;
	double tau;
	int iterationNumber;
	int particlesNumber;
	int allParticlesNumber;
	double U0;
	double density0;
	double B0;
	double upstreamR;
	double downstreamR;
	double deltaR;
	double temperature;
	double gridParameter;
	double epsilonR;
	double simulationTime;
	double deltaT;
	int A;
	int Z;

	int simulationType;

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

	double mass;

	int rgridNumber;

	double R0;
	double time;

	double maxSoundSpeed;

	Particle* Particles;
	SpaceBin** bins;

	double* grid;
	double* pointDensity;
	double* middleDensity;
	double* pointVelocity;
	double* middleMomentum;
	double* pointPressure;
	double* middleEnergy;

	double momentumAtPoint(int i);
	double energyAtPoint(int i);

	Simulation();
	~Simulation();

	void initializeProfile();
	void simulate();
	void evaluateHydrodynamic();
	void TracPen(double* u, double* flux, double cs, double leftFlux, double rightFlux);
	void updateValues();
	double minmod(double a, double b);
	void updateMaxSoundSpeed();

	void updateParameters();
};

#endif