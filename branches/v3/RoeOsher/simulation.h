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
	int prevShockWavePoint;
	bool shockWaveMoved;
	double shockWaveSpeed;

	double shockWaveT;

	double R0;
	double myTime;

	double maxSoundSpeed;

	double mass;
	double totalMomentum;
	double totalEnergy;
	double totalKineticEnergy;
	double totalTermalEnergy;
	double totalMagneticEnergy;
	double totalParticleEnergy;
	double totalParticles;
	double injectedParticles;

	double minP;
	double maxP;
	double* pgrid;
	double* logPgrid;
	double* kgrid;
	double* logKgrid;
	double deltaLogP;
	double deltaLogK;

	double* grid;
	double* middleGrid;
	double* deltaR;
	double* middleDeltaR;
	double* tempGrid;
	double* pointDensity;
	double* middleDensity;
	double* pointVelocity;
	double* middleVelocity;
	double* pointEnthalpy;
	double* middlePressure;
	double* cosmicRayPressure;
	double* pointSoundSpeed;
	double* tempDensity;
	double* tempMomentum;
	double* tempEnergy;

	double* dFlux;
	double* mFlux;
	double* eFlux;
	double** dFluxPlus;
	double** mFluxPlus;
	double** eFluxPlus;
	double** dFluxMinus;
	double** mFluxMinus;
	double** eFluxMinus;

	double* tempU;

	double** diffusionCoef;

	double** distributionFunction;
	double** tempDistributionFunction;

	double** magneticField;
	double** tempMagneticField;
	double** largeScaleField;
	double** growth_rate;
	double** crflux;
	double* integratedFlux;
	double* maxRate;

	double* magneticInductionSum;

	double momentum(int i);
	double energy(int i);
	double kineticEnergy(int i);
	double termalEnergy(int i);
	double temperatureIn(int i);
	double soundSpeed(int i);
	double volume(int i);
	double densityFlux(int i);
	void evaluateFluxes();
	void updateFluxes();
	void updateFluxes(double* flux, double** fluxPlus, double** fluxMinus);

	void updateDiffusionCoef();

	Simulation();
	~Simulation();

	void initializeProfile();
	void simulate();
	void evaluateHydrodynamic();
	void solveDiscontinious();
	void CheckNegativeDensity();
	void TracPen(double* u, double* flux, double cs);
	void updateFlux(double* flux);
	bool CheckShockWave(double& u, double p1, double p2,double u1, double u2, double rho1, double rho2);

	void evaluateCR();
	void solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x, double* alpha, double* beta);
	double injection(int i);
	void evaluateCosmicRayPressure();

	void evaluateField();
	void evaluateCRFlux();
	void growthRate();

	double minmod(double a, double b);
	double superbee(double a, double b);
	void updateMaxSoundSpeed();
	void updateShockWavePoint();
	void updateParameters();
	void updateTimeStep();

	void updateAll();

	void updateGrid();
	//std::list<GridZone*> createZones(int* type, double* gradientU, int& smallGradientZoneCount, int& bigGradientZoneCount);
    //void putPointsIntoZones(std::list<GridZone*>& zones, int pointsCount, int smallGradientZoneCount, int bigGradientZoneCount);
	//void convertZonesToGrid(std::list<GridZone*>& zones);
	//void addPoints(GridZone* zone, int& i);
	void redistributeValues();
	
};

#endif