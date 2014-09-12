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

	std::vector<MatrixElement>**** maxwellEquationMatrix;
	double**** maxwellEquationRightPart;

	Vector3d*** electricFlux;
	double*** electricDensity;
	Matrix3d*** dielectricTensor;
	Matrix3d*** pressureTensor;

	Vector3d*** Efield;
	Vector3d*** Bfield;

	Vector3d*** newEfield;
	Vector3d*** newBfield;

	Vector3d*** tempEfield;
	Vector3d*** tempBfield;

	std::vector<Particle*> particles;

	std::vector<Particle*>*** particlesInEbin;
	std::vector<Particle*>*** particlesInBbin;

	Matrix3d Kronecker;
	double LeviCivita[3][3][3];

	FILE* traectoryFile;
	FILE* distributionFile;
	FILE* EfieldFile;
	FILE* BfieldFile;
	FILE* Xfile;
	FILE* Yfile;
	FILE* Zfile;

	Simulation();
	Simulation(double xn, double yn, double zn, double xsizev, double ysizev, double zsizev, double temp, double rho, double Ex, double Ey, double Ez, double Bx, double By, double Bz, int maxIterations, double maxTimeV, int particlesPerBinV);
	~Simulation();

	void initialize();
	void createArrays();
	void createFiles();
	void simulate();

	void updateDeltaT();
	void createParticles();
	Particle* createParticle(int i, int j, int k, double weight, ParticleTypes type);

	Vector3d correlationTempEfield(Particle* particle);
	Vector3d correlationBfield(Particle* particle);
	Vector3d correlationTempEfield(Particle& particle);
	Vector3d correlationBfield(Particle& particle);
	Vector3d correlationEfield(Particle* particle);
	Vector3d correlationEfield(Particle& particle);
	
	Vector3d correlationFieldWithBbin(Particle& particle, int i, int j, int k);
	Vector3d correlationFieldWithEbin(Particle& particle, int i, int j, int k);
	Vector3d correlationFieldWithTempEbin(Particle& particle, int i, int j, int k);
	double correlationWithBbin(Particle& particle, int i, int j, int k);
	double correlationWithEbin(Particle& particle, int i, int , int k);
	double correlationBspline(const double& x, const double&  dx, const double& leftx, const double& rightx);

	Matrix3d Simulation::evaluateAlphaRotationTensor(double beta, Vector3d velocity, Vector3d EField, Vector3d BField); //see Noguchi

	void moveParticles();
	void correctParticlePosition(Particle* particle);
	void moveParticle(Particle* particle);
	void moveParticleNewtonIteration(Particle* particle, double* const oldCoordinates, double* const tempCoordinates, double* const newCoordinates);
	void evaluateParticlesRotationTensor();

	void evaluateFields();
	void checkMaxwellEquationMatrix();
	void createPerfectConductaryBoundaryCondition(int j, int k);
	void createInternalEquationX(int i, int j, int k);
	void createInternalEquationY(int i, int j, int k);
	void createInternalEquationZ(int i, int j, int k);
	void createInternalEquation(int i, int j, int k);
	void evaluateMaxwellEquationMatrix();
	void evaluateMagneticField();

	void generalizedMinimalResidualMethod();
	double***** arnoldiIterations(double** outHessenbergMatrix, int n, double***** prevBasis, double** prevHessenbergMatrix);
	double**** multiplySpecialMatrixVector(double**** vector);
	double evaluateError(double** hessenbergMatrix, double* vector, double beta, int n);
	double scalarMultiplyLargeVectors(double**** a, double**** b);

	double volume(int i, int j, int k);

	void collectParticlesIntoBins();
	void pushParticleIntoEbin(Particle* particle, int i, int j, int k);
	void pushParticleIntoBbin(Particle* particle, int i, int j, int k);
	bool particleCrossBbin(Particle& particle, int i, int j, int k);
	bool particleCrossEbin(Particle& particle, int i, int j, int k);
	void checkParticleInBox(Particle& particle);

	void updateParameters();
	void updateFields();
	double evaluateDivFlux(int i, int j, int k);
	Vector3d evaluateRotB(int i, int j, int k);
	Vector3d evaluateRotE(int i, int j, int k);
	Vector3d evaluateDivPressureTensor(int i, int j, int k);
	Vector3d evaluateGradDensity(int i, int j, int k);

	Vector3d getBfield(int i, int j, int k);
	Vector3d getTempEfield(int i, int j, int k);
	Vector3d getEfield(int i, int j, int k);
	Matrix3d getPressureTensor(int i, int j, int k);
	double getDensity(int i, int j, int k);
};

#endif