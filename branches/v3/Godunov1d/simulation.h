#ifndef SIMULATION_H
#define SIMULATION_H

#include <list>

class Simulation{
public:
	int iterationNumber;
	int particlesNumber;
	double U0;
	double density0;
	double B0;
	double upstreamR;
	double downstreamR;
	double temperature;
	double maxTime;
	double epsilonR;
	double deltaT;
	double deltaR0;
	double initialEnergy;
	int A;
	int Z;
	int currentIteration;

	int simulationType;
	bool tracPen;

	int rgridNumber;
	int shockWavePoint;
	bool shockWaveMoved;

	double R0;
	double myTime;

	double maxSoundSpeed;

	double mass;
	double totalMomentum;
	double totalEnergy;
	double totalKineticEnergy;
	double totalTermalEnergy;
	double totalParticles;
	double injectedParticles;

	double minP;
	double maxP;
	double* pgrid;
	double* logPgrid;
	double deltaLogP;

	double* grid;
	double* middleGrid;
	double* deltaR;
	double* middleDeltaR;
	double* tempGrid;
	double* pointDensity;
	double* middleDensity;
	double* pointVelocity;
	double* middleVelocity;
	double* pointPressure;
	double* middlePressure;
	double* cosmicRayPressure;

	double* distrFunDerivative;
	double* distrFunDerivative2;

	double* tempU;

	double** distributionFunction;
	double** tempDistributionFunction;

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

	double Simulation::diffusionCoef(int i, double p);

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
	void solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x, double* alpha, double* beta);
	double injection();
	void evaluateCosmicRayPressure();
	void changeDistrFunction();

	double minmod(double a, double b);
	void updateMaxSoundSpeed();
	void updateShockWavePoint();
	void updateParameters();

	void updateGrid();
	//std::list<GridZone*> createZones(int* type, double* gradientU, int& smallGradientZoneCount, int& bigGradientZoneCount);
    //void putPointsIntoZones(std::list<GridZone*>& zones, int pointsCount, int smallGradientZoneCount, int bigGradientZoneCount);
	//void convertZonesToGrid(std::list<GridZone*>& zones);
	//void addPoints(GridZone* zone, int& i);
	void redistributeValues();
	
};

#endif