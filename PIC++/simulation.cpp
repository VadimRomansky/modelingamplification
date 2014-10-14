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
	plasma_period2 = plasma_period;
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

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			delete[] EfieldX[i][j];
		}
		delete[] EfieldX[i];
	}
	delete[] EfieldX;

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			delete[] EfieldY[i][j];
		}
		delete[] EfieldY[i];
	}
	delete[] EfieldY;

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			delete[] EfieldZ[i][j];
		}
		delete[] EfieldZ[i];
	}
	delete[] EfieldZ;

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			delete[] BfieldX[i][j];
		}
		delete[] BfieldX[i];
	}
	delete[] BfieldX;

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			delete[] BfieldY[i][j];
		}
		delete[] BfieldY[i];
	}
	delete[] BfieldY;

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			delete[] BfieldZ[i][j];
		}
		delete[] BfieldZ[i];
	}
	delete[] BfieldZ;

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			delete[] electronConcentration[i][j];
			delete[] protonConcentration[i][j];
			delete[] chargeDensity[i][j];
		}
		delete[] electronConcentration[i];
		delete[] protonConcentration[i];
		delete[] chargeDensity[i];
	}

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

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				EfieldX[i][j][k] = E0.x;
				electricFluxX[i][j][k] = 0;
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				EfieldY[i][j][k] = E0.y;
				electricFluxY[i][j][k] = 0;
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber; ++k) {
				EfieldZ[i][j][k] = E0.z;
				electricFluxZ[i][j][k] = 0;
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				BfieldX[i][j][k] = B0.x;
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber; ++k) {
				BfieldY[i][j][k] = B0.y;
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				BfieldZ[i][j][k] = B0.z;
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				electronConcentration[i][j][k] = 0;
				protonConcentration[i][j][k] = 0;
				chargeDensity[i][j][k] = 0;
			}
		}
	}

	//electricDensity[xnumber/2][ynumber/2][znumber/2] = 1;
}

void Simulation::initializeSimpleElectroMagneticWave() {
	E0 = Vector3d(0, 0, 0);
	B0 = Vector3d(0, 0, 0);

	double kw = 2 * pi / xsize;
	double E = 1E-5;

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				EfieldX[i][j][k] = 0;
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				EfieldY[i][j][k] = E*sin(kw*xgrid[i]);
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber; ++k) {
				EfieldZ[i][j][k] = 0;
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				BfieldX[i][j][k] = 0;
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber; ++k) {
				BfieldY[i][j][k] = 0;
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				BfieldZ[i][j][k] = 0;
			}
		}
	}
}

void Simulation::createArrays() {
	xgrid = new double[xnumber + 1];
	ygrid = new double[ynumber + 1];
	zgrid = new double[znumber + 1];

	middleXgrid = new double[xnumber];
	middleYgrid = new double[ynumber];
	middleZgrid = new double[znumber];

	EfieldX = new double**[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		EfieldX[i] = new double*[ynumber + 1];
		for(int j = 0; j < ynumber + 1; ++j) {
			EfieldX[i][j] = new double[znumber + 1];
		}
	}

	EfieldY = new double**[xnumber+1];
	for(int i = 0; i < xnumber + 1; ++i) {
		EfieldY[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			EfieldY[i][j] = new double[znumber+1];
		}
	}

	EfieldZ = new double**[xnumber + 1];
	for(int i = 0; i < xnumber + 1; ++i) {
		EfieldZ[i] = new double*[ynumber + 1];
		for(int j = 0; j < ynumber + 1; ++j) {
			EfieldZ[i][j] = new double[znumber];
		}
	}

	electricFluxX = new double**[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		electricFluxX[i] = new double*[ynumber + 1];
		for(int j = 0; j < ynumber + 1; ++j) {
			electricFluxX[i][j] = new double[znumber + 1];
		}
	}

	electricFluxY = new double**[xnumber+1];
	for(int i = 0; i < xnumber + 1; ++i) {
		electricFluxY[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			electricFluxY[i][j] = new double[znumber+1];
		}
	}

	electricFluxZ = new double**[xnumber + 1];
	for(int i = 0; i < xnumber + 1; ++i) {
		electricFluxZ[i] = new double*[ynumber + 1];
		for(int j = 0; j < ynumber + 1; ++j) {
			electricFluxZ[i][j] = new double[znumber];
		}
	}

	BfieldX = new double**[xnumber + 1];
	for(int i = 0; i < xnumber + 1; ++i) {
		BfieldX[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			BfieldX[i][j] = new double[znumber];
		}
	}

	BfieldY = new double**[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		BfieldY[i] = new double*[ynumber + 1];
		for(int j = 0; j < ynumber + 1; ++j) {
			BfieldY[i][j] = new double[znumber];
		}
	}

	BfieldZ = new double**[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		BfieldZ[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			BfieldZ[i][j] = new double[znumber + 1];
		}
	}

	electronConcentration = new double**[xnumber];
	protonConcentration = new double**[xnumber];
	chargeDensity = new double**[xnumber];

	for(int i = 0; i < xnumber; ++i) {
		electronConcentration[i] = new double*[ynumber];
		protonConcentration[i] = new double*[ynumber];
		chargeDensity[i] = new double*[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			electronConcentration[i][j] = new double[znumber];
			protonConcentration[i][j] = new double[znumber];
			chargeDensity[i][j] = new double[znumber];
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
	//createParticles();
	updateDensityParameters();
	updateEnergy();

	//updateDeltaT();
	//deltaT = 0;

	while (time*plasma_period < maxTime && currentIteration < maxIteration) {
		printf("iteration number %d time = %15.10g\n", currentIteration, time * plasma_period);

		if (currentIteration % writeParameter == 0) {
			output();
		}
		updateDeltaT();

		moveParticles();
		updateElectroMagneticParameters();
		evaluateFields();
		
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
	outputFields(EfieldFile, BfieldFile, EfieldX, EfieldY, EfieldZ, BfieldX, BfieldY, BfieldZ, xnumber, ynumber, znumber, plasma_period, gyroradius);
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
	outputConcentrations(densityFile, electronConcentration, protonConcentration, chargeDensity, xnumber, ynumber, znumber, plasma_period, gyroradius);
	fclose(densityFile);

	divergenceErrorFile = fopen("./output/divergence_error.dat", "a");
	outputDivergenceError(divergenceErrorFile, this);
	fclose(divergenceErrorFile);
}

void Simulation::updateDeltaT() {
	printf("updating time step\n");
	double delta = min3(deltaX, deltaY, deltaZ);
	deltaT = 0.1 * delta / speed_of_light_normalized;
	double B = B0.norm();
	double E = E0.norm();
	/*for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				if(Bfield[i][j][k].norm() > B) {
					B = Bfield[i][j][k].norm();
				}
			}
		}
	}*/
	/*for(int i = 0; i < xnumber+1; ++i) {
		for(int j = 0; j < ynumber+1; ++j) {
			for(int k = 0; k < znumber+1; ++k) {
				if(Efield[i][j][k].norm() > E) {
					E = Efield[i][j][k].norm();
				}
			}
		}
	}*/
	if(B > 0){
		deltaT = min2(deltaT, 0.005 * massElectron * speed_of_light_normalized / (electron_charge_normalized * B));
	}
	if(E > 0) {
		deltaT = min2(deltaT, 0.005*massElectron * speed_of_light_normalized/(electron_charge_normalized*E));
	}
	//deltaT = 0.005 * massElectron * speed_of_light_normalized / (electron_charge_normalized * B0.norm());
	//deltaT = min2(deltaT, 0.02);
	deltaT = min2(deltaT, plasma_period2/10);
	//deltaT = min2(deltaT, 1E-1);
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

	double dx = deltaX / 4;
	double dy = deltaY / 4;
	double dz = deltaZ / 4;

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
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k){
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
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