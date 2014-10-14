#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "stdlib.h"
#include "stdio.h"
#include "vector"
#include "matrix3d.h"
#include "matrixElement.h"
#include "particle.h"

class Simulation{
public:
	int xnumber;
	int ynumber;
	int znumber;
	int particlesNumber;
	int particlesPerBin;

	double density;
	double temperature;

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
	Simulation(double xn, double yn, double zn, double xsizev, double ysizev, double zsizev, double temp, double rho, double Ex, double Ey, double Ez, double Bx, double By, double Bz, int maxIterations, double maxTimeV, int particlesPerBinV);
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

	double correlationBspline(const double& x, const double&  dx, const double& leftx, const double& rightx);

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

	void updateEnergy();

	double volume(int i, int j, int k);

	void updateElectroMagneticParameters();
	void updateDensityParameters();
};

#endif