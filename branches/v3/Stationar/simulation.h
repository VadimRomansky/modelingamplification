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
	double initialEnergy;
	int currentIteration;

	int rgridNumber;
	int shockWavePoint;
	int prevShockWavePoint;

	double R0;

	double maxSoundSpeed;

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
	double* vscattering;

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
	double* cosmicRayConcentration;

	double** magneticField;
	double** tempMagneticField;
	double** largeScaleField;
	double* magneticEnergy;
	double* tempMagneticEnergy;
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

	void evaluateCR();
	void solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x, double* alpha, double* beta);
	double injection(int i);
	void evaluateCosmicRayPressure();

	void evaluateField();
	void evaluateCRFlux();
	void growthRate();

	double minmod(double a, double b);
	double superbee(double a, double b);

	void updateAll();
};

#endif