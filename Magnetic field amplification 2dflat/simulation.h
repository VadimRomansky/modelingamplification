#ifndef SIMULATION_H
#define SIMULATION_H

#include "stdafx.h"
#include "particle.h"
#include "SpaceBin.h"
#include "constants.h"
#include "vector3d.h"

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
	double X1;
	double X2;
	double deltaX;
	double Y1;
	double Y2;
	double deltaY;
	double temperature;

	double epsilonR;
	int A;
	int Z;

	bool smallAngleScattering;

	double energy;
	double theorEnergy;

	double momentumX;
	double theorMomentumX;

	double momentumY;
	double theorMomentumY;

	double particlesWeight;

	int xgridNumber;
	int ygridNumber;
	int shockWaveIndex;
	double shockWavePoint;
	int currentShockWavePoint;

	int freeTimeEvaluationType;

	Particle* Particles;
	SpaceBin*** bins;

	std::list <Particle*> startPDF;
	std::vector<Particle*> introducedParticles;
	std::vector<Particle*> escapedParticles;

	Simulation();
	~Simulation();

	void initializeProfile();
	void simulate();
	void resetDetectors();
	std::vector<Particle*> getParticles();
	void resetProfile();
	void collectAverageVelocity();
	void removeEscapedParticles();
	void updateEnergy();
	void updateShockWavePoint();
	void updateCosmicRayBoundMomentum(bool write);

	//void findShockWavePoint();

private:
	static const int particleMultiply = 20;
};



#endif