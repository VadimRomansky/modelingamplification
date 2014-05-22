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

	double velocityCD;
	double velocityL1;
	double velocityL2;
	double velocityR1;
	double velocityR2;

	double rCD;
	double rL1;
	double rL2;
	double rR1;
	double rR2;

	double p1;
	double rho1;
	double u1;
	double p2;
	double rho2;
	double u2;
	double p3;
	double rho3;
	double u3;
	double p4;
	double rho4;
	double u4;
	double p5;
	double rho5;
	double u5;
	double p6;
	double rho6;
	double u6;

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

	Simulation();
	~Simulation();

	void initializeProfile();
	void simulate();
	void evaluateHydrodynamic();
	void solveDiscontinious();
	void CheckNegativeDensity();
	void TracPen(double* u, double* flux, double cs);
	void updateFlux(double* flux);

	void riemanSolve();
	void riemanMove();
	double pressureFunction(double p, double p1, double rho1);
	double pressureFunctionDerivative(double p, double p1, double rho1);
	void successiveApproximationPressure(double& p, double& u, double& R1, double& R2, double& alpha1, double& alpha2, double p1, double p2, double u1, double u2, double rho1, double rho2);
	double firstApproximationPressure(double rho1, double rho2, double u1, double u2, double p1, double p2);

	double minmod(double a, double b);
	double superbee(double a, double b);
	void updateMaxSoundSpeed();
	void updateShockWavePoint();
	void updateParameters();
	void updateTimeStep();

	void updateAll();
	
};

#endif