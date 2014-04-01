#ifndef SIMULATION_H
#define SIMULATION_H

#include <list>

class Simulation{
public:
	int iterationNumber;
	double U0;
	double density0;
	double B0;
	double upstreamR;
	double downstreamR;
	double temperature;
	double maxTime;
	double deltaT;
	double deltaR0;
	double deltaLogP;

	int currentIteration;

	int rgridNumber;
	int shockWavePoint;

	double R0;
	double myTime;

	double totalParticles;
	double injectedParticles;

	double minP;
	double maxP;
	double* pgrid;
	double* logPgrid;

	double* grid;
	double* middleGrid;
	double* deltaR;
	double* middleDeltaR;


	double** distributionFunction;
	double** tempDistributionFunction;

	double velocity(double x);
	double density(double x);
	double volume(int i);

	double diffusionCoef(double x, double p);

	Simulation();
	~Simulation();

	void simulate();
	void initializeProfile();

	void evaluateCR();
	void solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x, double* alpha, double* beta);
	double injection();

	void updateParameters();
	
};

#endif