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

	int simulationType;
	bool tracPen;

	int rgridNumber;
	int shockWavePoint;

	double R0;
	double time;

	double maxSoundSpeed;

	double mass;
	double totalMomentum;
	double totalEnergy;
	double totalKineticEnergy;
	double totalTermalEnergy;

	double minP;
	double maxP;
	double* pgrid;

	double* grid;
	double* middleGrid;
	double* pointDensity;
	double* middleDensity;
	double* pointVelocity;
	double* middleVelocity;
	double* pointPressure;
	double* middlePressure;

	double** distributionFunction;
	double** tempDistribution;
	double** distributionDerivative;

	double momentum(int i);
	double energy(int i);
	double kineticEnergy(int i);
	double termalEnergy(int i);
	double temperatureIn(int i);
	double soundSpeed(int i);
	double volume(int i);

	double densityFlux(int i);
	double momentumConvectiveFlux(int i);
	double energyFlux(int i);

	double diffussionCoef(int i, int j);

	Simulation();
	~Simulation();

	void initializeProfile();
	void simulate();
	void evaluateHydrodynamic();
	void solveDiscontinious();
	double pressureFunction(double p, double p1, double rho1);
	double pressureFunctionDerivative(double p, double p1, double rho1);
	double pressureFunctionDerivative2(double p, double p1, double rho1);
	void successiveApproximationPressure(double& p, double& u, double& R1, double& R2, double& alpha1, double& alpha2, double p1, double p2, double u1, double u2, double rho1, double rho2);
	double firstApproximationPressure(double rho1, double rho2, double u1, double u2, double p1, double p2);
	void CheckNegativeDensity();
	void TracPen(double* u, double* flux, double cs);

	void evaluateCR();
	void solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x);
	double injection();

	double minmod(double a, double b);
	void updateMaxSoundSpeed();
	void updateShockWavePoint();

	void updateParameters();
};

#endif