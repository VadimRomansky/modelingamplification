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
	double upstreamR;
	double downstreamR;
	double deltaR;
	double deltaTheta;
	double deltaPhi;
	double temperature;
	double gridParameter;
	double alpha;
	double beta;
	double gamma;
	double delta;
	double epsilonR;
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

	double momentumX;
	double theorMomentumX;
	double particlesWeight;

	int rgridNumber;
	int shockWaveIndex;
	double shockWavePoint;
	int currentShockWavePoint;

	int simulationType;
	int freeTimeEvaluationType;

	Particle* Particles;
	SpaceBin** bins;
	SpaceBin* zeroBin;

	double zeroBinScale;
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
	void introduceNewParticles();
	void updateMaxMinK();
	void updateMaxMinP(double& minP, double& maxP);
	void evaluateMagneticField(double* startField,double* endField,double deltax,int gridNumber,double gradU, double density,double U,int binNumber);
	SpaceBin* getStartBin(double theta, double phi);
	void detectFromTo(int fromIndexR, int fromIndexTheta, int fromIndexPhi, int toIndexR, int toIndexTheta, int toIndexPhi, const Particle& particle);
	Particle* getAnyParticle();
	void resetProfile();
	void collectAverageVelocity();
	void sortParticlesIntoBins();
	void removeEscapedParticles();
	void updateEnergy();
	void updateShockWavePoint();
	void updateCosmicRayBoundMomentum(bool write);
	void smoothProfile();
	void smoothProfile(std::list<SpaceBin*> bins);

	void findShockWavePoint();
	void saveEnergyAndMomentum();

private:
	static const int particleMultiply = 20;
};



#endif