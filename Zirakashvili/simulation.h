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

	double mass;
	double downstreamMass1;
	double downstreamMass2;
	double upstreamMass1;
	double upstreamMass2;

	int rgridNumber;

	double R0;
	double time;

	double forwardShockWaveR;
	double oldForwardShockWaveR;
	double reverseShockWaveR;
	double oldReverseShockWaveR;
	double contactDiscontR;
	double oldContactDiscontR;

	double forwardV;
	double reverseV;
	double contactDiscontV;

	double maxSoundSpeed;

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
	void moveShockWaves();
	void updateDownstreamValues();
	void TracPen(double* u, double* flux, double cs, double deltaXi, double leftFlux, double rightFlux);
	double minmod(double a, double b);
	void updateBinsVolume();
	void updateMaxSoundSpeed();

	void updateParameters();
};

#endif