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

	int simulationType;
	bool tracPen;

	int rgridNumber;

	int shockWavePoint;
	int prevShockWavePoint;
	bool shockWaveMoved;
	double shockWaveSpeed;

	double shockWaveT;

	double R0;
	double myTime;

	int currentIteration;

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

	double deltaLogK;

	double* grid;
	double* gridsquare;
	double* volumeFactor;
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
	double* gridVelocity;
	double* tempDensity;
	double* tempMomentum;
	double* tempEnergy;

	double* kgrid;
	double* logKgrid;

	double** diffusionCoef;
	double* tempU;

	double** distributionFunction;
	double** tempDistributionFunction;

	double** magneticField;
	double** tempMagneticField;
	double** largeScaleField;
	double** growth_rate;
	double** crflux;

	double* magneticInductionSum;

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
	void TracPenRadial(double* u, double* flux, double cs);

	void evaluateCR();
	void solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x, double* alpha, double* beta);
	double injection();

	double minmod(double a, double b);
	void updateMaxSoundSpeed();
	void updateShockWavePoint();
	void updateParameters();
	void updateTimeStep();
	void evaluateCosmicRayPressure();

	void updateDiffusionCoef();
	void evaluateField();
	void evaluateCRFlux();
	void growthRate();

	void updateGrid();

	void updateAll();
	//std::list<GridZone*> createZones(int* type, double* gradientU, int& smallGradientZoneCount, int& bigGradientZoneCount);
    //void putPointsIntoZones(std::list<GridZone*>& zones, int pointsCount, int smallGradientZoneCount, int bigGradientZoneCount);
	//void convertZonesToGrid(std::list<GridZone*>& zones);
	//void addPoints(GridZone* zone, int& i);
	void redistributeValues();
	
};

#endif