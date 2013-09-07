#ifndef SIMULATION_H
#define SIMULATION_H

#include "stdafx.h"
#include "SpaceBin.h"
#include "particle.h"

class Simulation{
public:
	int iterationNumber;
	int particlesNumber;
	double U0;
	double density0;
	double B0;
	double upstreamR;
	double downstreamR;
	double deltaR;
	double temperature;
	double epsilonR;
	double deltaT;
	int A;
	int Z;

	bool cycleBound;
	int simulationType;

	double mass;

	int rgridNumber;

	double R0;
	double time;

	double maxSoundSpeed;

	double* grid;
	double* pointDensityLeft;
	double* pointDensityRight;
	double* middleDensity;
	double* pointVelocity;
	double* middleVelocity;
	double* pointPressure;
	double* middlePressure;

	double momentum(int i);
	double energy(int i);
	double temperatureIn(int i);
	double soundSpeed(int i);

	Simulation();
	~Simulation();

	void initializeProfile();
	void simulate();
	void evaluateHydrodynamic();
	void solveDiscontinious();
	void evaluatePressureAndVelocity();
	void TracPen(double* u, double* flux, double cs, double leftFlux, double rightFlux);
	void updateValues();
	double minmod(double a, double b);
	void updateMaxSoundSpeed();

	void updateParameters();
};

#endif