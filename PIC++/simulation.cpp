#include "stdlib.h"
#include "stdio.h"

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "output.h"
#include "constants.h"
#include "matrix3d.h"
#include "random.h"

Simulation::Simulation() {
	debugMode = true;

	currentIteration = 0;
	time = 0;
	particlesNumber = 0;

	theta = 0.5;

	particleEnergy = 0;
	electricFieldEnergy = 0;
	magneticFieldEnergy = 0;

	momentum = Vector3d(0, 0, 0);

	//read input!
	xnumber = 10;
	ynumber = 10;
	znumber = 10;

	xsize = 1E15;
	ysize = 1E15;
	zsize = 1E15;

	temperature = 0.1E14;
	density = 1.6E-24;

	maxIteration = 1E7;
	maxTime = 1E7;

	particlesPerBin = 10;

	B0 = Vector3d(1E-5, 0, 0);
	E0 = Vector3d(1E-15, 0, 0);

	double concentration = density / massProton;

	plasma_period = sqrt(massElectron / (4 * pi * concentration * sqr(electron_charge))) / (2 * pi);

	double thermal_momentum;
	if(kBoltzman*temperature > massElectron*speed_of_light*speed_of_light) {
		thermal_momentum = kBoltzman*temperature/speed_of_light;
	} else {
		thermal_momentum = sqrt(2*massElectron*kBoltzman*temperature);
	}
	gyroradius = thermal_momentum*speed_of_light / (electron_charge * B0.norm());

	kBoltzman_normalized = kBoltzman * gyroradius * gyroradius / (plasma_period * plasma_period);
	speed_of_light_normalized = speed_of_light * plasma_period / gyroradius;
	speed_of_light_normalized_sqr = speed_of_light_normalized * speed_of_light_normalized;
	electron_charge_normalized = electron_charge * plasma_period / sqrt(cube(gyroradius));

	E0 = E0 / (plasma_period * gyroradius);
	B0 = B0 / (plasma_period * gyroradius);

	density = density * cube(gyroradius);

	Kronecker = Matrix3d(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0);

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				LeviCivita[i][j][k] = 0;
			}
		}
	}
	LeviCivita[0][1][2] = 1.0;
	LeviCivita[0][2][1] = -1.0;
	LeviCivita[1][0][2] = -1.0;
	LeviCivita[1][2][0] = 1.0;
	LeviCivita[2][0][1] = 1.0;
	LeviCivita[2][1][0] = -1.0;
}

Simulation::Simulation(double xn, double yn, double zn, double xsizev, double ysizev, double zsizev, double temp, double rho, double Ex, double Ey, double Ez, double Bx, double By, double Bz, int maxIterations, double maxTimeV, int particlesPerBinV) {
	debugMode = true;

	currentIteration = 0;
	time = 0;
	particlesNumber = 0;

	particleEnergy = 0;
	electricFieldEnergy = 0;
	magneticFieldEnergy = 0;

	momentum = Vector3d(0, 0, 0);


	theta = 0.5;

	//read input!
	xnumber = xn;
	ynumber = yn;
	znumber = zn;

	xsize = xsizev;
	ysize = ysizev;
	zsize = zsizev;

	temperature = temp;
	density = rho;

	maxIteration = maxIterations;
	maxTime = maxTimeV;

	particlesPerBin = particlesPerBinV;

	B0 = Vector3d(Bx, By, Bz);
	E0 = Vector3d(Ex, Ey, Ez);

	double concentration = density / massProton;

	plasma_period = sqrt(massElectron / (4 * pi * concentration * sqr(electron_charge))) / (2 * pi);
	double thermal_momentum;
	if(kBoltzman*temperature > massElectron*speed_of_light*speed_of_light) {
		thermal_momentum = kBoltzman*temperature/speed_of_light;
	} else {
		thermal_momentum = sqrt(2*massElectron*kBoltzman*temperature);
	}
	gyroradius = thermal_momentum*speed_of_light / (electron_charge * B0.norm());

	plasma_period = 1.0;
	gyroradius = 1.0;

	kBoltzman_normalized = kBoltzman * plasma_period * plasma_period / (gyroradius * gyroradius);
	speed_of_light_normalized = speed_of_light * plasma_period / gyroradius;
	speed_of_light_normalized_sqr = speed_of_light_normalized * speed_of_light_normalized;
	electron_charge_normalized = electron_charge * plasma_period / sqrt(cube(gyroradius));

	E0 = E0 * (plasma_period * gyroradius);
	B0 = B0 * (plasma_period * gyroradius);

	density = density * cube(gyroradius);

	xsize /= gyroradius;
	ysize /= gyroradius;
	zsize /= gyroradius;
	printf("xsize/gyroradius = %lf\n", xsize);

	Kronecker = Matrix3d(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0);

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				LeviCivita[i][j][k] = 0;
			}
		}
	}
	LeviCivita[0][1][2] = 1.0;
	LeviCivita[0][2][1] = -1.0;
	LeviCivita[1][0][2] = -1.0;
	LeviCivita[1][2][0] = 1.0;
	LeviCivita[2][0][1] = 1.0;
	LeviCivita[2][1][0] = -1.0;
}

Simulation::~Simulation() {
	for (int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		delete particle;
	}
	particles.clear();

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				particlesInBbin[i][j][k].clear();
			}
			delete[] Bfield[i][j];
			delete[] newBfield[i][j];
			//delete[] tempBfield[i][j];
			delete[] electricDensity[i][j];
			delete[] pressureTensor[i][j];

			delete[] electronConcentration[i][j];
			delete[] protonConcentration[i][j];
			delete[] chargeDensity[i][j];
		}
		delete[] Bfield[i];
		delete[] newBfield[i];
		//delete[] tempBfield[i];
		delete[] electricDensity[i];
		delete[] pressureTensor[i];

		delete[] electronConcentration[i];
		delete[] protonConcentration[i];
		delete[] chargeDensity[i];

	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				particlesInEbin[i][j][k].clear();
			}
			delete[] Efield[i][j];
			delete[] newEfield[i][j];
			delete[] tempEfield[i][j];
			delete[] electricFlux[i][j];
			delete[] dielectricTensor[i][j];
		}
		delete[] Efield[i];
		delete[] newEfield[i];
		delete[] tempEfield[i];
		delete[] electricFlux[i];
		delete[] dielectricTensor[i];
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					maxwellEquationMatrix[i][j][k][l].clear();
				}
				delete[] maxwellEquationMatrix[i][j][k];
				delete[] maxwellEquationRightPart[i][j][k];

			}
			delete[] maxwellEquationMatrix[i][j];
			delete[] maxwellEquationRightPart[i][j];
		}
		delete[] maxwellEquationMatrix[i];
		delete[] maxwellEquationRightPart[i];
	}
	delete[] maxwellEquationMatrix;
	delete[] maxwellEquationRightPart;

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				delete[] divergenceCleaningPotential[i][j][k];
				delete[] divergenceCleanUpMatrix[i][j][k];
				delete[] divergenceCleanUpRightPart[i][j][k];
			}
			delete[] divergenceCleaningPotential[i][j];
			delete[] divergenceCleanUpMatrix[i][j];
			delete[] divergenceCleanUpRightPart[i][j];
		}
		delete[] divergenceCleaningPotential[i];
		delete[] divergenceCleanUpMatrix[i];
		delete[] divergenceCleanUpRightPart[i];
	}
	delete[] divergenceCleaningPotential;
	delete[] divergenceCleanUpMatrix;
	delete[] divergenceCleanUpRightPart;

	delete[] Efield;
	delete[] newEfield;
	delete[] tempEfield;
	delete[] Bfield;
	delete[] newBfield;
	//delete[] tempBfield;

	delete[] electricDensity;
	delete[] electricFlux;
	delete[] dielectricTensor;
	delete[] pressureTensor;

	delete[] electronConcentration;
	delete[] protonConcentration;
	delete[] chargeDensity;

	delete[] xgrid;
	delete[] ygrid;
	delete[] zgrid;
	delete[] middleXgrid;
	delete[] middleYgrid;
	delete[] middleZgrid;
}

void Simulation::initialize() {

	deltaX = xsize / (xnumber);
	deltaY = ysize / (ynumber);
	deltaZ = zsize / (znumber);

	for (int i = 0; i <= xnumber; ++i) {
		xgrid[i] = i * deltaX;
	}

	for (int i = 0; i <= ynumber; ++i) {
		ygrid[i] = i * deltaY;
	}

	for (int i = 0; i <= znumber; ++i) {
		zgrid[i] = i * deltaZ;
	}

	for (int i = 0; i < xnumber; ++i) {
		middleXgrid[i] = (xgrid[i] + xgrid[i + 1]) / 2;
	}

	for (int i = 0; i < ynumber; ++i) {
		middleYgrid[i] = (ygrid[i] + ygrid[i + 1]) / 2;
	}

	for (int i = 0; i < znumber; ++i) {
		middleZgrid[i] = (zgrid[i] + zgrid[i + 1]) / 2;
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k] = E0;
				newEfield[i][j][k] = Efield[i][j][k];
				tempEfield[i][j][k] = Efield[i][j][k];

			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Bfield[i][j][k] = B0;
				Bfield[i][j][k].x = B0.x;
				newBfield[i][j][k] = Bfield[i][j][k];
				//tempBfield[i][j][k] = Bfield[i][j][k];
				electronConcentration[i][j][k] = 0;
				protonConcentration[i][j][k] = 0;
				chargeDensity[i][j][k] = 0;
				electricDensity[i][j][k] = 0;
			}
		}
	}
}

void Simulation::initializeSimpleElectroMagneticWave() {
	E0 = Vector3d(0, 0, 0);
	B0 = Vector3d(0, 0, 0);
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Bfield[i][j][k] = Vector3d(0, 0, 0);
				newBfield[i][j][k] = Bfield[i][j][k];
			}
		}
	}
	double kw = 2 * pi / xsize;
	double E = 1E-5;

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k].x = 0;
				Efield[i][j][k].y = E * sin(kw * xgrid[i]);
				Efield[i][j][k].z = 0;
				tempEfield[i][j][k] = Efield[i][j][k];
				newEfield[i][j][k] = tempEfield[i][j][k];
			}
		}
	}

	double t = 2 * pi / (kw * speed_of_light_normalized);
}

void Simulation::createArrays() {
	xgrid = new double[xnumber + 1];
	ygrid = new double[ynumber + 1];
	zgrid = new double[znumber + 1];

	middleXgrid = new double[xnumber];
	middleYgrid = new double[ynumber];
	middleZgrid = new double[znumber];

	Efield = new Vector3d**[xnumber + 1];
	newEfield = new Vector3d**[xnumber + 1];
	tempEfield = new Vector3d**[xnumber + 1];
	Bfield = new Vector3d**[xnumber];
	newBfield = new Vector3d**[xnumber];
	//tempBfield = new Vector3d**[xnumber];

	electronConcentration = new double**[xnumber];
	protonConcentration = new double**[xnumber];
	chargeDensity = new double**[xnumber];

	electricFlux = new Vector3d**[xnumber + 1];
	electricDensity = new double**[xnumber];
	dielectricTensor = new Matrix3d**[xnumber + 1];
	pressureTensor = new Matrix3d**[xnumber];

	divergenceCleaningPotential = new double***[xnumber+1];

	particlesInEbin = new std::vector<Particle*>**[xnumber + 1];
	particlesInBbin = new std::vector<Particle*>**[xnumber];

	for (int i = 0; i < xnumber; ++i) {
		Bfield[i] = new Vector3d*[ynumber];
		newBfield[i] = new Vector3d*[ynumber];
		//tempBfield[i] = new Vector3d*[ynumber];
		electricDensity[i] = new double*[ynumber];
		pressureTensor[i] = new Matrix3d*[ynumber];
		particlesInBbin[i] = new std::vector<Particle*>*[ynumber];

		electronConcentration[i] = new double*[ynumber];
		protonConcentration[i] = new double*[ynumber];
		chargeDensity[i] = new double*[ynumber];

		for (int j = 0; j < ynumber; ++j) {
			Bfield[i][j] = new Vector3d[znumber];
			newBfield[i][j] = new Vector3d[znumber];
			//tempBfield[i][j] = new Vector3d[znumber];
			electricDensity[i][j] = new double[znumber];
			pressureTensor[i][j] = new Matrix3d[znumber];
			particlesInBbin[i][j] = new std::vector<Particle*>[znumber];

			electronConcentration[i][j] = new double[znumber];
			protonConcentration[i][j] = new double[znumber];
			chargeDensity[i][j] = new double[znumber];
		}
	}

	for (int i = 0; i < xnumber + 1; ++i) {
		Efield[i] = new Vector3d*[ynumber + 1];
		newEfield[i] = new Vector3d*[ynumber + 1];
		tempEfield[i] = new Vector3d*[ynumber + 1];
		electricFlux[i] = new Vector3d*[ynumber + 1];
		dielectricTensor[i] = new Matrix3d*[ynumber + 1];
		particlesInEbin[i] = new std::vector<Particle*>*[ynumber + 1];

		for (int j = 0; j < ynumber + 1; ++j) {
			Efield[i][j] = new Vector3d[znumber + 1];
			newEfield[i][j] = new Vector3d[znumber + 1];
			tempEfield[i][j] = new Vector3d[znumber + 1];
			electricFlux[i][j] = new Vector3d[znumber + 1];
			dielectricTensor[i][j] = new Matrix3d[znumber + 1];
			particlesInEbin[i][j] = new std::vector<Particle*>[znumber + 1];
		}
	}

	maxwellEquationMatrix = new std::vector<MatrixElement>***[xnumber];
	maxwellEquationRightPart = new double***[xnumber];
	for (int i = 0; i < xnumber; ++i) {
		maxwellEquationMatrix[i] = new std::vector<MatrixElement>**[ynumber];
		maxwellEquationRightPart[i] = new double**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			maxwellEquationMatrix[i][j] = new std::vector<MatrixElement>*[znumber];
			maxwellEquationRightPart[i][j] = new double*[znumber];
			for (int k = 0; k < znumber; ++k) {
				maxwellEquationMatrix[i][j][k] = new std::vector<MatrixElement>[3];
				maxwellEquationRightPart[i][j][k] = new double[3];
			}
		}
	}

	divergenceCleanUpMatrix = new std::vector<MatrixElement>***[xnumber+1];
	divergenceCleanUpRightPart = new double***[xnumber+1];

	for(int i = 0; i < xnumber + 1; ++i) {
		divergenceCleaningPotential[i] = new double**[ynumber];
		divergenceCleanUpMatrix[i] = new std::vector<MatrixElement>**[ynumber];
		divergenceCleanUpRightPart[i] = new double**[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			divergenceCleaningPotential[i][j] = new double*[znumber];
			divergenceCleanUpMatrix[i][j] = new std::vector<MatrixElement>*[znumber];
			divergenceCleanUpRightPart[i][j] = new double*[znumber];
			for(int k = 0; k < znumber; ++k) {
				divergenceCleaningPotential[i][j][k] = new double[1];
				divergenceCleanUpMatrix[i][j][k] = new std::vector<MatrixElement>[1];
				divergenceCleanUpRightPart[i][j][k] = new double[1];
			}
		}
	}
}

void Simulation::createFiles() {
	protonTraectoryFile = fopen("./output/traectory_proton.dat", "w");
	fclose(protonTraectoryFile);
	electronTraectoryFile = fopen("./output/traectory_electron.dat", "w");
	fclose(electronTraectoryFile);
	distributionFile = fopen("./output/distribution_protons.dat", "w");
	fclose(distributionFile);
	EfieldFile = fopen("./output/Efield.dat", "w");
	fclose(EfieldFile);
	BfieldFile = fopen("./output/Bfield.dat", "w");
	fclose(BfieldFile);
	Xfile = fopen("./output/Xfile.dat", "w");
	fclose(Xfile);
	Yfile = fopen("./output/Yfile.dat", "w");
	fclose(Yfile);
	Zfile = fopen("./output/Zfile.dat", "w");
	fclose(Zfile);
	generalFile = fopen("./output/general.dat", "w");
	fclose(generalFile);
	densityFile = fopen("./output/concentrations.dat", "w");
	fclose(densityFile);
	divergenceErrorFile = fopen("./output/divergence_error.dat", "w");
	fclose(divergenceErrorFile);
}

void Simulation::simulate() {
	createArrays();
	initialize();
	initializeSimpleElectroMagneticWave();
	createFiles();
	createParticles();
	collectParticlesIntoBins();
	updateDensityParameters();
	//cleanupDivergence();
	updateEnergy();

	//updateDeltaT();
	//deltaT = 0;

	while (time*plasma_period < maxTime && currentIteration < maxIteration) {
		printf("iteration number %d time = %15.10g\n", currentIteration, time * plasma_period);

		if (currentIteration % writeParameter == 0) {
			output();
		}
		updateDeltaT();
		evaluateParticlesRotationTensor();
		evaluateFields();
		moveParticles();
		updateFields();
		//cleanupDivergence();
		updateDensityParameters();
		updateEnergy();

		time += deltaT;
		currentIteration++;

	}
}

void Simulation::output() {
	printf("outputing\n");
	if (particles.size() > 0) {
		distributionFile = fopen("./output/distribution_protons.dat", "a");
		outputDistribution(distributionFile, particles, ParticleTypes::PROTON);
		fclose(distributionFile);
		protonTraectoryFile = fopen("./output/traectory_proton.dat", "a");
		outputTraectory(protonTraectoryFile, getFirstProton(), time);
		fclose(protonTraectoryFile);
		electronTraectoryFile = fopen("./output/traectory_electron.dat", "a");
		outputTraectory(electronTraectoryFile, getFirstElectron(), time);
		fclose(electronTraectoryFile);
	}
	EfieldFile = fopen("./output/Efield.dat", "a");
	BfieldFile = fopen("./output/Bfield.dat", "a");
	outputFields(EfieldFile, BfieldFile, Efield, Bfield, xnumber, ynumber, znumber, plasma_period, gyroradius);
	fclose(EfieldFile);
	fclose(BfieldFile);

	Xfile = fopen("./output/Xfile.dat", "w");
	outputGrid(Xfile, xgrid, xnumber);
	fclose(Xfile);

	Yfile = fopen("./output/Yfile.dat", "w");
	outputGrid(Yfile, ygrid, ynumber);
	fclose(Yfile);

	Zfile = fopen("./output/Zfile.dat", "w");
	outputGrid(Zfile, zgrid, znumber);
	fclose(Zfile);

	generalFile = fopen("./output/general.dat", "a");
	outputGeneral(generalFile, this);
	fclose(generalFile);

	densityFile = fopen("./output/concentrations.dat", "a");
	outputConcentrations(densityFile, electronConcentration, protonConcentration, chargeDensity, electricDensity, xnumber, ynumber, znumber, plasma_period, gyroradius);
	fclose(densityFile);

	divergenceErrorFile = fopen("./output/divergence_error.dat", "a");
	outputDivergenceError(divergenceErrorFile, this);
	fclose(divergenceErrorFile);
}


Matrix3d Simulation::evaluateAlphaRotationTensor(double beta, Vector3d velocity, Vector3d EField, Vector3d BField) {
	Matrix3d result;

	beta = beta / speed_of_light_normalized;

	double gamma_factor = 1 / sqrt(1 - (velocity.x * velocity.x + velocity.y * velocity.y + velocity.z * velocity.z) / speed_of_light_normalized_sqr);
	double G = ((beta * (EField.scalarMult(velocity)) / speed_of_light_normalized_sqr) + gamma_factor);
	beta = beta / G;
	double denominator = G * (1 + beta * beta * BField.scalarMult(BField));

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			result.matrix[i][j] = Kronecker.matrix[i][j] + beta * beta * BField[i] * BField[j];
			for (int k = 0; k < 3; ++k) {
				for (int l = 0; l < 3; ++l) {
					result.matrix[i][j] -= beta * LeviCivita[j][k][l] * Kronecker.matrix[i][l] * BField[k];
				}
			}

			result.matrix[i][j] /= denominator;
		}
	}

	return result;
}

void Simulation::updateDeltaT() {
	printf("updating time step\n");
	double delta = min3(deltaX, deltaY, deltaZ);
	deltaT = 0.1 * delta / speed_of_light_normalized;
	double B = B0.norm();
	double E = E0.norm();
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				if(Bfield[i][j][k].norm() > B) {
					B = Bfield[i][j][k].norm();
				}
			}
		}
	}
	for(int i = 0; i < xnumber+1; ++i) {
		for(int j = 0; j < ynumber+1; ++j) {
			for(int k = 0; k < znumber+1; ++k) {
				if(Efield[i][j][k].norm() > E) {
					E = Efield[i][j][k].norm();
				}
			}
		}
	}
	if(B > 0){
		deltaT = min2(deltaT, 0.005 * massElectron * speed_of_light_normalized / (electron_charge_normalized * B));
	}
	if(E > 0) {
		deltaT = min2(deltaT, 0.005*massElectron * speed_of_light_normalized/(electron_charge_normalized*E));
	}
	//deltaT = 0.005 * massElectron * speed_of_light_normalized / (electron_charge_normalized * B0.norm());
	//deltaT = min2(deltaT, 0.02);
	//deltaT = min2(deltaT, plasma_period/10);
	deltaT = min2(deltaT, 1E-1);
}

void Simulation::createParticles() {
	printf("creating particles\n");
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				double weight = (density / (massProton * particlesPerBin)) * volume(i, j, k);
				//double weight = (1.0 / particlesPerBin) * volume(i, j, k);
				Vector3d coordinates;
				for (int l = 0; l < 2 * particlesPerBin; ++l) {
					ParticleTypes type;
					if (l % 2 == 0) {
						type = ParticleTypes::PROTON;
					} else {
						type = ParticleTypes::ELECTRON;
					}
					Particle* particle = createParticle(i, j, k, weight, type);
					if (l % 2 == 0) {
						coordinates = particle->coordinates;
					} else {
						particle->coordinates= coordinates;
					}
					particles.push_back(particle);
					particlesNumber++;
					if(particlesNumber % 1000 == 0){
						printf("create particle number %d\n", particlesNumber);
					}
				}
			}
		}
	}
}

Particle* Simulation::getFirstProton() {
	for(int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		if(particle->type == ParticleTypes::PROTON) {
			return particle;
		}
	}
	return NULL;
}

Particle* Simulation::getFirstElectron() {
	for(int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		if(particle->type == ParticleTypes::ELECTRON) {
			return particle;
		}
	}
	return NULL;
}

Particle* Simulation::createParticle(int i, int j, int k, double weight, ParticleTypes type) {
	double charge = 0;
	double mass = 0;

	switch (type) {
	case ParticleTypes::PROTON:
		mass = massProton;
		charge = electron_charge_normalized;
		break;
	case ParticleTypes::ELECTRON:
		mass = massElectron;
		charge = -electron_charge_normalized;
		break;
	}

	double x = xgrid[i] + deltaX * uniformDistribution();
	double y = ygrid[j] + deltaY * uniformDistribution();
	double z = zgrid[k] + deltaZ * uniformDistribution();

	double dx = deltaX / 100;
	double dy = deltaY / 100;
	double dz = deltaZ / 100;

	double energy = mass * speed_of_light_normalized_sqr;

	double thetaParamter = kBoltzman_normalized * temperature / (mass * speed_of_light_normalized_sqr);

	if (thetaParamter < 0.01) {
		energy = mass * speed_of_light_normalized_sqr + maxwellDistribution(temperature, kBoltzman_normalized);
	} else {
		energy = maxwellJuttnerDistribution(temperature, mass, speed_of_light_normalized, kBoltzman_normalized);
	}

	double p = sqrt(energy * energy - sqr(mass * speed_of_light_normalized_sqr)) / speed_of_light_normalized;

	double pz = p * (2 * uniformDistribution() - 1);
	double phi = 2 * pi * uniformDistribution();
	double pnormal = sqrt(p * p - pz * pz);
	double px = pnormal * cos(phi);
	double py = pnormal * sin(phi);

	Particle* particle = new Particle(mass, charge, weight, type, x, y, z, px, py, pz, dx, dy, dz);

	return particle;
}

double Simulation::volume(int i, int j, int k) {
	return deltaX * deltaY * deltaZ;
}

void Simulation::checkParticleInBox(Particle& particle) {
	if (particle.coordinates.x < 0) {
		printf("particle.x < 0\n");
		exit(0);
	}
	if (particle.coordinates.x > xgrid[xnumber]) {
		printf("particle.x > xsize\n");
		exit(0);
	}
	if (particle.coordinates.y < 0) {
		printf("particle.y < 0\n");
		exit(0);
	}
	if (particle.coordinates.y > ygrid[ynumber]) {
		printf("particle.y > ysize\n");
		exit(0);
	}
	if (particle.coordinates.z < 0) {
		printf("particle.z < 0\n");
		exit(0);
	}
	if (particle.coordinates.z > zgrid[znumber]) {
		printf("particle.z > zsize\n");
		exit(0);
	}
}

void Simulation::updateElectroMagneticParameters() {
	printf("updating flux, density snd dielectric tensor\n");
	collectParticlesIntoBins();
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				electricFlux[i][j][k] = Vector3d(0, 0, 0);
				dielectricTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				for (int pcount = 0; pcount < particlesInEbin[i][j][k].size(); ++pcount) {
					Particle* particle = particlesInEbin[i][j][k][pcount];
					double correlation = correlationWithEbin(*particle, i, j, k) / volume(i, j, k);
					if(i == 0) {
						correlation = correlation * 2; //because smaller volume
					}
					Vector3d velocity = particle->velocity(speed_of_light_normalized);
					double gamma = particle->gammaFactor(speed_of_light_normalized);
					Vector3d rotatedVelocity = particle->rotationTensor * velocity * gamma;

					electricFlux[i][j][k] += rotatedVelocity * particle->charge * particle->weight * correlation;
					dielectricTensor[i][j][k] = dielectricTensor[i][j][k] - particle->rotationTensor * (theta * deltaT * deltaT * 2 * pi * particle->charge * particle->charge * correlation / particle->mass);
				}
			}
		}
	}

	//for i = 0
	for(int j = 0; j < ynumber + 1; ++j) {
		for(int k = 0; k < znumber + 1; ++k) {
			double rho = 0;
			for(int pcount = 0; pcount < particlesInEbin[0][j][k].size(); ++pcount) {
				Particle* particle = particlesInEbin[0][j][k][pcount];
				Particle tempParticle = Particle(*particle);

				double correlation = 2*correlationWithEbin(tempParticle, -1, j, k)/volume(0, j, k);
				tempParticle.momentum.x = - tempParticle.momentum.x;
				double beta = tempParticle.charge*deltaT/tempParticle.mass;
				Vector3d velocity = tempParticle.velocity(speed_of_light_normalized);

				Vector3d oldE = correlationEfield(particle);
				Vector3d oldB = correlationBfield(particle);

				tempParticle.rotationTensor = evaluateAlphaRotationTensor(beta, velocity, oldE, oldB);

				double gamma = tempParticle.gammaFactor(speed_of_light_normalized);
				Vector3d rotatedVelocity = tempParticle.rotationTensor * velocity * gamma;

				electricFlux[0][j][k] += rotatedVelocity * tempParticle.charge * tempParticle.weight * correlation;
				dielectricTensor[0][j][k] = dielectricTensor[0][j][k] - tempParticle.rotationTensor * (theta * deltaT * deltaT * 2 * pi * tempParticle.charge * tempParticle.charge * correlation / tempParticle.mass);

			}
		}
	}

	//for periodic conditions we must summ sides parameters
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			electricFlux[i][j][0] = electricFlux[i][j][0] + electricFlux[i][j][znumber];
			electricFlux[i][j][znumber] = electricFlux[i][j][0];

			dielectricTensor[i][j][0] = dielectricTensor[i][j][0] + dielectricTensor[i][j][znumber];
			dielectricTensor[i][j][znumber] = dielectricTensor[i][j][0];
		}

		for (int k = 0; k < znumber + 1; ++k) {
			electricFlux[i][0][k] = electricFlux[i][0][k] + electricFlux[i][ynumber][k];
			electricFlux[i][ynumber][k] = electricFlux[i][0][k];

			dielectricTensor[i][0][k] = dielectricTensor[i][0][k] + dielectricTensor[i][ynumber][k];
			dielectricTensor[i][ynumber][k] = dielectricTensor[i][0][k];
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				electricDensity[i][j][k] = 0;
				pressureTensor[i][j][k] = Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0);
				for (int pcount = 0; pcount < particlesInBbin[i][j][k].size(); ++pcount) {
					Particle* particle = particlesInBbin[i][j][k][pcount];
					double correlation = correlationWithBbin(*particle, i, j, k) / volume(i, j, k);

					double gamma = particle->gammaFactor(speed_of_light_normalized);
					Vector3d velocity = particle->velocity(speed_of_light_normalized);
					Vector3d rotatedVelocity = particle->rotationTensor * velocity * gamma;

					electricDensity[i][j][k] += particle->weight * particle->charge * correlation;

					pressureTensor[i][j][k] = rotatedVelocity.tensorMult(rotatedVelocity) * particle->weight * particle->charge * correlation;
				}
			}
		}
	}

	//for i = 0
	for(int j = 0; j < ynumber; ++j) {
		for(int k = 0; k < znumber; ++k) {
			double rho = 0;
			for(int pcount = 0; pcount < particlesInBbin[0][j][k].size(); ++pcount) {
				Particle* particle = particlesInBbin[0][j][k][pcount];
				Particle tempParticle = Particle(*particle);

				double correlation = correlationWithBbin(tempParticle, -1, j, k)/volume(0, j, k);
				tempParticle.momentum.x = - tempParticle.momentum.x;
				double beta = tempParticle.charge*deltaT/tempParticle.mass;
				Vector3d velocity = tempParticle.velocity(speed_of_light_normalized);

				Vector3d oldE = correlationEfield(particle);
				Vector3d oldB = correlationBfield(particle);

				tempParticle.rotationTensor = evaluateAlphaRotationTensor(beta, velocity, oldE, oldB);

				double gamma = tempParticle.gammaFactor(speed_of_light_normalized);
				Vector3d rotatedVelocity = tempParticle.rotationTensor * velocity * gamma;

				electricDensity[0][j][k] += tempParticle.weight*tempParticle.charge*correlation;
				pressureTensor[0][j][k] += rotatedVelocity.tensorMult(rotatedVelocity)*tempParticle.weight*tempParticle.charge*correlation;

			}
		}
	}

	for (int i = 0; i <= xnumber; ++i) {
		for (int j = 0; j <= ynumber; ++ j) {
			for (int k = 0; k <= znumber; ++k) {
				Vector3d divPressureTensor = evaluateDivPressureTensor(i, j, k);
				electricFlux[i][j][k] = electricFlux[i][j][k] - divPressureTensor * deltaT / 2;
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				double divJ = evaluateDivFlux(i, j, k);

				electricDensity[i][j][k] -= deltaT * theta * divJ;
			}
		}
	}
}

void Simulation::updateDensityParameters() {
	double full_density = 0;
	double full_p_concentration = 0;
	double full_e_concentration = 0;
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				electronConcentration[i][j][k] = 0;
				protonConcentration[i][j][k] = 0;
				chargeDensity[i][j][k] = 0;

				for(int pcount = 0; pcount < particlesInBbin[i][j][k].size(); ++pcount) {
					Particle* particle = particlesInBbin[i][j][k][pcount];

					double correlation = correlationWithBbin(*particle, i, j, k)/volume(i, j, k);
					if(i == 0) {
						correlation += correlationWithBbin(*particle, i-1, j, k)/volume(i-1, j, k);
					}

					chargeDensity[i][j][k] += correlation*particle->charge*particle->weight;
					if(particle->type == ParticleTypes::ELECTRON) {
						electronConcentration[i][j][k] += correlation*particle->weight;
					} else if (particle->type == ParticleTypes::PROTON) {
						protonConcentration[i][j][k] += correlation*particle->weight;
					}
				}
				if(i < xnumber - 1  && i > 0){
					full_density += chargeDensity[i][j][k]*volume(i, j, k);
					full_p_concentration += protonConcentration[i][j][k]*volume(i, j, k);
					full_e_concentration += electronConcentration[i][j][k]*volume(i, j, k);
					electronConcentration[i][j][k] /= cube(gyroradius);
					protonConcentration[i][j][k] /= cube(gyroradius);
					chargeDensity[i][j][k] /= (sqrt(cube(gyroradius))*plasma_period);
				}
			}
		}
	}

	full_density/= ((xsize-2*deltaX)*ysize*zsize);
	full_p_concentration /= ((xsize-2*deltaX)*ysize*zsize*cube(gyroradius));
	full_e_concentration /= ((xsize-2*deltaX)*ysize*zsize*cube(gyroradius));
}

void Simulation::updateEnergy() {
	particleEnergy = 0;
	electricFieldEnergy = 0;
	magneticFieldEnergy = 0;
	energy = 0;

	momentum = Vector3d(0, 0, 0);
	for(int i = 0; i < xnumber+1; ++i) {
		for(int j = 0; j < ynumber+1; ++j){
			for(int k = 0; k < znumber+1; ++k){
				double factor = 1;
				if(i == 0 || i == xnumber) {
					factor = factor/2;
				}
				if(j == 0 || j == ynumber) {
					factor = factor/2;
				}
				if(k == 0 || k == znumber) {
					factor = factor/2;
				}
				Vector3d E = Efield[i][j][k];
				electricFieldEnergy += E.scalarMult(E)*volume(i, j, k)*factor/(8*pi);
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k){
				Vector3d B = Bfield[i][j][k];
				magneticFieldEnergy += B.scalarMult(B)*volume(i, j, k)/(8*pi);
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				Vector3d E = (Efield[i][j][k]  + Efield[i][j+1][k] + Efield[i][j][k+1] + Efield[i][j+1][k+1] + Efield[i+1][j][k] + Efield[i+1][j+1][k] + Efield[i+1][j][k+1] + Efield[i+1][j+1][k+1])/8;
				Vector3d B = Bfield[i][j][k];

				momentum += (E.vectorMult(B)/(4*pi*speed_of_light_normalized))*volume(i, j, k);
			}
		}
	}

	for(int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];

		particleEnergy += particle->energy(speed_of_light_normalized)*particle->weight;
		momentum += particle->momentum*particle->weight;
	}

	particleEnergy *= sqr(gyroradius/plasma_period);
	electricFieldEnergy *= sqr(gyroradius/plasma_period);
	magneticFieldEnergy *= sqr(gyroradius/plasma_period);
	momentum = momentum * gyroradius/plasma_period;


	energy = particleEnergy + electricFieldEnergy + magneticFieldEnergy;
}

Vector3d Simulation::getBfield(int i, int j, int k) {
	if (i < 0) {
		return Vector3d(0.0, 0.0, 0.0);
	} else if (i >= xnumber) {
		return B0;
	}

	if (j < 0) {
		j = ynumber - 1;
	} else if (j >= ynumber) {
		j = 0;
	}

	if (k < 0) {
		k = znumber - 1;
	} else if (k >= znumber) {
		k = 0;
	}

	return Bfield[i][j][k];
}

Vector3d Simulation::getTempEfield(int i, int j, int k) {
	if (i < 0) {
		//return Vector3d(0.0, 0.0, 0.0);
		i = 0;
	} else if (i > xnumber) {
		return E0;
	}

	if (j < 0) {
		j = j + ynumber;
	} else if (j > ynumber) {
		j = j - ynumber;
	}

	if (k < 0) {
		k = k + znumber;
	} else if (k > znumber) {
		k = k - znumber;
	}

	return tempEfield[i][j][k];
}

Vector3d Simulation::getEfield(int i, int j, int k) {
	if (i < 0) {
		return Vector3d(0.0, 0.0, 0.0);
	} else if (i > xnumber) {
		return E0;
	}

	if (j < 0) {
		j = ynumber - 1;
	} else if (j > ynumber) {
		j = 1;
	}

	if (k < 0) {
		k = znumber - 1;
	} else if (k > znumber) {
		k = 1;
	}

	return Efield[i][j][k];
}

Matrix3d Simulation::getPressureTensor(int i, int j, int k) {
	if (i < 0) {
		return Matrix3d(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	} else if (i >= xnumber) {
		return Matrix3d(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	}

	if (j < 0) {
		j = ynumber - 1;
	} else if (j >= ynumber) {
		j = 0;
	}

	if (k < 0) {
		k = znumber - 1;
	} else if (k >= znumber) {
		k = 0;
	}

	return pressureTensor[i][j][k];
}

double Simulation::getDensity(int i, int j, int k) {
	if (i < 0) return 0;
	if (i >= xnumber) return 0;

	if (j < 0) {
		j = ynumber - 1;
	}
	if (j >= ynumber) {
		j = 0;
	}

	if (k < 0) {
		k = znumber - 1;
	}
	if (k >= znumber) {
		k = 0;
	}

	return electricDensity[i][j][k];
}