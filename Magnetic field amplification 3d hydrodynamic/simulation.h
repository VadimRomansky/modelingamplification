#ifndef SIMULATION_H
#define SIMULATION_H

#include "stdafx.h"
//#include "particle.h"
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
	//int cosmicRayNumber;
	double partOfCosmicRay;
	double minK;
	double maxK;

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
	void resetDetectors();
	//std::list <Particle> getParticleUniformDistribution(double pmax,int pnumber,int thetanumber); 
	std::list<Particle> getParticleGaussDistribution(int number);
	std::vector<Particle*> getParticles();
	double maxwell(double p, double temperature,double mass);
	void updateMagneticField();
	vector3d gradientSpeed(int i, int j, int k);
	double cascadingDerivativeW(double w,double k,double rho);
	double derivativeFieldK(double* field,int j);
	double instability(double w,int i, int j);
	double cascadingDerivativeK(double w,double k,double rho);
	double dissipation(double w, int i, double k);
	void updateMaxMinK();
	void updateMaxMinP(double& minP, double& maxP);
	void updatePressureSpectralDensity();
	void evaluateMagneticField(double* startField,double* endField,double deltax,int gridNumber,double gradU, double density,double U,int binNumber);
	SpaceBin* getStartBin(double theta, double phi);
	void detectFromTo(int fromIndexR, int fromIndexTheta, int fromIndexPhi, int toIndexR, int toIndexTheta, int toIndexPhi, const Particle& particle);
	void multiplyParticleWeights(double v);
	Particle* getAnyParticle();
	void resetProfile();
	void collectAverageVelocity();
	void sortParticlesIntoBins();
	void removeEscapedParticles();
	void updateEnergy();
	void updateShockWavePoint();
	void updateCosmicRayBoundMomentum();
	void evaluateHydrodynamic(double* newDensity, double* newMomentum, double* newPressure);


	double densityFlux(int i);
	double momentumFlux(int i);
	double pressureFlux(int i);
	double momentumPressureFlux(int i);
	double volumeFlux(int i);

	double vanleer(double a, double b);
	double findMaxVelocity();

	void tvd(double* value, double* flux, double maxVelocity);
	void convectionTVD(double* value, double* flux);
private:
	static const int particleMultiply = 20;
};



#endif