#ifndef XBIN_H
#define XBIN_H

#include "stdafx.h"
#include "constants.h"
#include "random.h"
//#include "matrix3d.h"
//#include "particle2.h"
#include <list>

class Particle;

class SpaceBin{
public:
	static const int defaultTimeRelation = 5;

	double x;
	double x1;
	double x2;

	double y;
	double y1;
	double y2;

	double Ux;
	double Uy;
	double density;
	double volume;
	double temperature;
	double B;
	double B0;
	double W;

	double initialMomentum;

	int numberX;
	int numberY;

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
	double averageVelocityX;
	double averageVelocityY;

	double crFlux;

	double centralMomentum;
	double momentumDifference;
	double energyDifference;

	int freeTimeEvaluationType;

	SpaceBin();
	SpaceBin(double x, double deltax, double y, double deltay, double ux, double uy, double rho, double temperature, double b, int i, int j, bool scattering, int freeTimeEvaluation);
	~SpaceBin();
	int* propagateParticle(Particle* particle ,double& time, double timeStep, const int xgridNumber, const int ygridNumber);
	//double getFreePath(Particle* particle);
	double getFreeTime(Particle* particle);
	void makeOneStep(Particle* particle, double deltat, double& time);
	bool isInBin(Particle* particle);
	bool isInThisOrNear(double r);
	static int* binByCoordinates(double x, double y, double x0, double deltax, const int xgridNumber, double y0, double deltay, const int ygridNumber); 
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
	double getU();

	void detectParticleR1(Particle* particle);
	void detectParticleR2(Particle* particle);

	void updateCosmicRayBoundMomentum(bool write);
	void updateTemperature(double* distribution, double deltap);

	void largeAngleScattering(Particle* particle, double& time, double timeStep);

	void evaluateU(double ux1, double ux2, double uy1, double uy2);
	double evaluateCrymskyIntegralX();
	double evaluateCrymskyIntegralY();

	//void operator=(Xbin bin);

};


#endif