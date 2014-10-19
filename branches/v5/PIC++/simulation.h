#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "stdlib.h"
#include "stdio.h"
#include "vector"
#include "matrix3d.h"
#include "matrixElement.h"
#include "particle.h"

enum BoundaryConditionTypes {SUPERCONDUCTERLEFT, PERIODIC};

class Simulation{
public:
	int xnumber;
	int ynumber;
	int znumber;
	int particlesNumber;
	int particlesPerBin;

	double density;
	double temperature;

	BoundaryConditionTypes boundaryConditionType;

	double plasma_period;
	double plasma_period2;
	double gyroradius;

	double speed_of_light_normalized;
	double speed_of_light_normalized_sqr;
	double kBoltzman_normalized;
	double electron_charge_normalized;

	double time;
	double maxTime;
	int currentIteration;
	int maxIteration;
	double xsize;
	double ysize;
	double zsize;

	double deltaT;

	double deltaX;
	double deltaY;
	double deltaZ;

	double theta;

	bool debugMode;

	double particleEnergy;
	double electricFieldEnergy;
	double magneticFieldEnergy;
	double energy;

	Vector3d momentum;

	double*** electronConcentration;
	double*** protonConcentration;
	double*** chargeDensity;


	Vector3d V0;

	Vector3d B0;
	Vector3d E0;
	Matrix3d pressureTensor0;

	double* xgrid;
	double* ygrid;
	double* zgrid;
	double* middleXgrid;
	double* middleYgrid;
	double* middleZgrid;

	double*** EfieldX;
	double*** EfieldY;
	double*** EfieldZ;

	double*** electricFluxX;
	double*** electricFluxY;
	double*** electricFluxZ;

	double*** BfieldX;
	double*** BfieldY;
	double*** BfieldZ;

	double*** divergenceCleaningPotential;
	std::vector<MatrixElement>*** divergenceCleanUpMatrix;
	double*** divergenceCleanUpRightPart;


	std::vector<Particle*> particles;

	Matrix3d Kronecker;
	double LeviCivita[3][3][3];

	FILE* protonTraectoryFile;
	FILE* electronTraectoryFile;
	FILE* distributionFile;
	FILE* EfieldFile;
	FILE* BfieldFile;
	FILE* Xfile;
	FILE* Yfile;
	FILE* Zfile;
	FILE* generalFile;
	FILE* densityFile;
	FILE* divergenceErrorFile;

	Simulation();
	Simulation(double xn, double yn, double zn, double xsizev, double ysizev, double zsizev, double temp, double rho, double Ex, double Ey, double Ez, double Bx, double By, double Bz, int maxIterations, double maxTimeV, int particlesPerBinV, BoundaryConditionTypes type);
	~Simulation();

	void initialize();
	void initializeSimpleElectroMagneticWave();
	void createArrays();
	void createFiles();
	void simulate();
	void output();

	void updateDeltaT();
	void createParticles();
	Particle* createParticle(int i, int j, int k, double weight, ParticleTypes type);
	Particle* getFirstProton();
	Particle* getFirstElectron();


	Vector3d correlationEfield(Particle* particle);
	Vector3d correlationBfield(Particle* particle);
	double correlationEfieldX(Particle* particle);
	double correlationEfieldY(Particle* particle);
	double correlationEfieldZ(Particle* particle);
	double correlationBfieldX(Particle* particle);
	double correlationBfieldY(Particle* particle);
	double correlationBfieldZ(Particle* particle);
	double correlationWithExBin(int i, int j, int k, Particle* particle);
	double correlationWithEyBin(int i, int j, int k, Particle* particle);
	double correlationWithEzBin(int i, int j, int k, Particle* particle);
	double correlationWithBxBin(int i, int j, int k, Particle* particle);
	double correlationWithByBin(int i, int j, int k, Particle* particle);
	double correlationWithBzBin(int i, int j, int k, Particle* particle);
	double correlationWithGridBin(int i, int j, int k, Particle* particle);
	double correlationWithShiftGridBin(int i, int j, int k, Particle*particle);

	double correlationX(int i, Particle* particle);
	double correlationShiftX(int i, Particle* particle);
	double correlationY(int j, Particle* particle);
	double correlationShiftY(int j, Particle* particle);
	double correlationZ(int k, Particle* particle);
	double correlationShiftZ(int k, Particle* particle);

	double correlationBspline(const double& x, const double&  dx, const double& leftx, const double& rightx);
	double getEx(int i, int j, int k);
	double getEy(int i, int j, int k);
	double getEz(int i, int j, int k);
	double getBx(int i, int j, int k);
	double getBy(int i, int j, int k);
	double getBz(int i, int j, int k);

	void moveParticles();
	void correctParticlePosition(Particle* particle);
	void moveParticle(Particle* particle);
	void checkParticleInBox(Particle& particle);

	void evaluateFields();
	void updateEfield(double dt);
	void updateEfieldX(int i, int j, int k, double dt);
	void updateEfieldY(int i, int j, int k, double dt);
	void updateEfieldZ(int i, int j, int k, double dt);
	void updateBfield(double dt);
	void updateBfieldX(int i, int j, int k, double dt);
	void updateBfieldY(int i, int j, int k, double dt);
	void updateBfieldZ(int i, int j, int k, double dt);
	void updateBoundaries();
	double evaluateDivE(int i, int j, int k);

	void cleanupDivergence();
	void updateFieldByPotential();
	double cleanUpRightPart(int i, int j, int k);
	void createDivergenceCleanupInternalEquation(int i, int j, int k);
	void createDivergenceCleanupRightEquation(int j, int k);
	void updateEnergy();

	double volume(int i, int j, int k);

	void resetElectroMagneticParameters();
	void updateElectroMagneticParameters();
	void updateElectricFluxX(Particle* particle);
	void updateElectricFluxY(Particle* particle);
	void updateElectricFluxZ(Particle* particle);
	void updateChargeDensity(Particle* particle);
	void updateBoundariesParameters();
	void updateElectroMagneticParameters(Particle* particle);
	void addElectricFluxX(int i, int j, int k, double flux);
	void addElectricFluxY(int i, int j, int k, double flux);
	void addElectricFluxZ(int i, int j, int k, double flux);
	void addChargeDensity(int i, int j, int k, double charge);
	void addConcentration(int i, int j, int k, double weight, ParticleTypes particle_type);

	double evaluateFullChargeDensity();
};

#endif