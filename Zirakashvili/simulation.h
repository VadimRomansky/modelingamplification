#ifndef SIMULATION_H
#define SIMULATION_H

#include "stdafx.h"
#include "SpaceBin.h"
#include "particle.h"

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
	double forwardShockWavePoint;
	double reverseShockWavePoint;
	double contactDiscontPoin;

	double xiMin;
	double xiMax;

	Particle* Particles;
	SpaceBin** upstreamBins1;
	SpaceBin** upstreamBins2;
	SpaceBin** downstreamBins1;
	SpaceBin** downstreamBins2;

    void solveUpstream1();
	void solveUpstream2();
	void solveDownstream1();
	void solveDownstream2();

	Simulation();
	~Simulation();

	void initializeProfile();
	void simulate();
};

#endif