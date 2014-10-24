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

	boundaryConditionType = SUPERCONDUCTERLEFT;

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

Simulation::Simulation(double xn, double yn, double zn, double xsizev, double ysizev, double zsizev, double temp, double rho, double Ex, double Ey, double Ez, double Bx, double By, double Bz, int maxIterations, double maxTimeV, int particlesPerBinV, BoundaryConditionTypes type) {
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

	boundaryConditionType = type;

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

	double lambda_debye = sqrt((kBoltzman*temperature)/(4*pi*electron_charge*electron_charge*concentration));

	double dx = xsize/xnumber;
	double dy = ysize/ynumber;
	double dz = zsize/znumber;

	if(particlesPerBin*cube(lambda_debye)/(dx*dy*dz) < 10)
	{
		printf("Number of macroparticles in debye sphere is to small\n");
	}

	if(dx > lambda_debye)
	{
		printf("dx > Debye length\n");
	}

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

	for (int i = 0; i < xnumber+1; ++i) {
		for (int j = 0; j < ynumber+1; ++j) {
			delete[] electronConcentration[i][j];
			delete[] protonConcentration[i][j];
			delete[] chargeDensity[i][j];
			delete[] velocity[i][j];
			delete[] divergenceCleaningPotential[i][j];
			delete[] divergenceCleanUpMatrix[i][j];
			delete[] divergenceCleanUpRightPart[i][j];
		}
		delete[] electronConcentration[i];
		delete[] protonConcentration[i];
		delete[] chargeDensity[i];
		delete[] velocity[i];
		delete[] divergenceCleaningPotential[i];
		delete[] divergenceCleanUpMatrix[i];
		delete[] divergenceCleanUpRightPart[i];
	}

	delete[] electronConcentration;
	delete[] protonConcentration;
	delete[] chargeDensity;
	delete[] velocity;
	delete[] divergenceCleaningPotential;
	delete[] divergenceCleanUpMatrix;
	delete[] divergenceCleanUpRightPart;

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

	for (int i = 0; i < xnumber+1; ++i) {
		for (int j = 0; j < ynumber+1; ++j) {
			for (int k = 0; k < znumber+1; ++k) {
				electronConcentration[i][j][k] = 0;
				protonConcentration[i][j][k] = 0;
				chargeDensity[i][j][k] = 0;
				velocity[i][j][k] = Vector3d(0, 0, 0);
				divergenceCleaningPotential[i][j][k] = 0;
				divergenceCleanUpRightPart[i][j][k] = 0;
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
	double omega = kw*speed_of_light_normalized;
	double t = 2*pi/omega;

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
				BfieldZ[i][j][k] = E*sin(kw*(xgrid[i]+deltaX/2));
			}
		}
	}
}

void Simulation::initializeAlfvenWave() {
	E0 = Vector3d(0, 0, 0);
	B0 = Vector3d(1E-3, 0, 0);

	double alfvenV = B0.norm()/sqrt(4*pi*density);
	if(alfvenV > speed_of_light_normalized) {
		printf("alfven velocity > c\n");
		exit(0);
	}

	double kw = 1 * 2 * pi / xsize;

	double omega = kw*alfvenV;
	double t = 2*pi/omega;
	double B = 1E-4;

	double qLog = 15;
	double nuElectronIon = 4*sqrt(2*pi)*qLog*power(electron_charge_normalized, 4)*(density/massProton)/(3*sqrt(massElectron)*sqrt(cube(temperature)));
	double collisionlessParameter = omega/nuElectronIon;
	double conductivity = (3*massProton*sqrt(massElectron)*sqrt(cube(kBoltzman_normalized*temperature)))/(sqr(massElectron)*4*sqrt(2*pi)*qLog*sqr(electron_charge_normalized));

	double alphaParameter = omega/conductivity;

	double magneticReynolds = conductivity*alfvenV/(kw*speed_of_light_normalized_sqr);

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
				EfieldY[i][j][k] =0;
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
				BfieldX[i][j][k] = B0.x;
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
				BfieldZ[i][j][k] = -B*cos(kw*(xgrid[i] + deltaX/2));
			}
		}
	}

	for(int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		Vector3d velocity = Vector3d(0, 0, 1)*(B/sqrt(4*pi*density))*cos(kw*particle->coordinates.x);
		double beta = velocity.norm()/speed_of_light_normalized;
		particle->addVelocity(velocity, speed_of_light_normalized);
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

	electronConcentration = new double**[xnumber+1];
	protonConcentration = new double**[xnumber+1];
	chargeDensity = new double**[xnumber+1];
	velocity = new Vector3d**[xnumber + 1];
	divergenceCleaningPotential = new double**[xnumber+1];
	divergenceCleanUpMatrix = new std::vector<MatrixElement>**[xnumber+1];
	divergenceCleanUpRightPart = new double**[xnumber + 1];

	for(int i = 0; i < xnumber+1; ++i) {
		electronConcentration[i] = new double*[ynumber+1];
		protonConcentration[i] = new double*[ynumber+1];
		chargeDensity[i] = new double*[ynumber+1];
		velocity[i] = new Vector3d*[ynumber + 1];
		divergenceCleaningPotential[i] = new double*[ynumber+1];
		divergenceCleanUpMatrix[i] = new std::vector<MatrixElement>*[ynumber+1];
		divergenceCleanUpRightPart[i] = new double*[ynumber + 1];
		for(int j = 0; j < ynumber+1; ++j) {
			electronConcentration[i][j] = new double[znumber+1];
			protonConcentration[i][j] = new double[znumber+1];
			chargeDensity[i][j] = new double[znumber+1];
			velocity[i][j] = new Vector3d[znumber + 1];
			divergenceCleaningPotential[i][j] = new double[znumber+1];
			divergenceCleanUpMatrix[i][j] = new std::vector<MatrixElement>[znumber+1];
			divergenceCleanUpRightPart[i][j] = new double[znumber + 1];
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
	velocityFile = fopen("./output/velocity.dat", "w");
	fclose(velocityFile);
	fluxFile = fopen("./output/flux.dat", "w");
	fclose(fluxFile);
	divergenceErrorFile = fopen("./output/divergence_error.dat", "w");
	fclose(divergenceErrorFile);
}

void Simulation::simulate() {
	createArrays();
	initialize();
	//initializeSimpleElectroMagneticWave();
	createFiles();
	createParticles();
	//initializeAlfvenWave();
	updateEnergy();
	updateElectroMagneticParameters();
	updateDeltaT();

	while (time*plasma_period < maxTime && currentIteration < maxIteration) {
		printf("iteration number %d time = %15.10g dt = %15.10g\n", currentIteration, time * plasma_period, deltaT * plasma_period);

		if (currentIteration % writeParameter == 0) {
			output();
		}
		updateDeltaT();

		moveParticles();
		updateElectroMagneticParameters();
		//evaluateFields();
		
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

	velocityFile = fopen("./output/velocity.dat", "a");
	outputVelocity(velocityFile, velocity, xnumber, ynumber, znumber, plasma_period, gyroradius);
	fclose(velocityFile);

	fluxFile = fopen("./output/flux.dat", "a");
	outputFlux(fluxFile, electricFluxX, electricFluxY, electricFluxZ, xnumber, ynumber, znumber, plasma_period, gyroradius);
	fclose(fluxFile);

	divergenceErrorFile = fopen("./output/divergence_error.dat", "a");
	outputDivergenceError(divergenceErrorFile, this);
	fclose(divergenceErrorFile);
}

void Simulation::updateDeltaT() {
	printf("updating time step\n");
	double delta = min3(deltaX, deltaY, deltaZ);
	deltaT = 0.3 * delta / speed_of_light_normalized;
	double B = B0.norm();
	double E = E0.norm();
	Particle* minElectron = electronMinMomentum();
	double minMomentum = minElectron->momentum.norm();
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				double Bfield = sqrt(BfieldX[i][j][k]*BfieldX[i][j][k] + BfieldY[i][j][k]*BfieldY[i][j][k] + BfieldZ[i][j][k]*BfieldZ[i][j][k]);
				if (Bfield > B)
				{
					B = Bfield;
				}
			}
		}
	}
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				double Efield = sqrt(EfieldX[i][j][k]*EfieldX[i][j][k] + EfieldY[i][j][k]*EfieldY[i][j][k] + EfieldZ[i][j][k]*EfieldZ[i][j][k]);
				if (Efield > E)
				{
					E = Efield;
				}
			}
		}
	}
	if(B > 0){
		deltaT = min2(deltaT, 0.02 * massElectron * speed_of_light_normalized / (electron_charge_normalized * B));
	}
	if(E > 0) {
		deltaT = min2(deltaT, 0.02* massElectron * speed_of_light_normalized/(electron_charge_normalized*E));
	}

	//if(E > 0) {
		//deltaT = min2(deltaT, 0.02*minMomentum/(electron_charge_normalized*E));
	//}
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
						type = PROTON;
					} else {
						type = ELECTRON;
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
		if(particle->type == PROTON) {
			return particle;
		}
	}
	return NULL;
}

Particle* Simulation::getFirstElectron() {
	for(int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		if(particle->type == ELECTRON) {
			return particle;
		}
	}
	return NULL;
}

Particle* Simulation::createParticle(int i, int j, int k, double weight, ParticleTypes type) {
	double charge = 0;
	double mass = 0;

	switch (type) {
	case PROTON:
		mass = massProton;
		charge = electron_charge_normalized;
		break;
	case ELECTRON:
		mass = massElectron;
		charge = -electron_charge_normalized;
		break;
	}

	double x = xgrid[i] + deltaX * uniformDistribution();
	double y = ygrid[j] + deltaY * uniformDistribution();
	double z = zgrid[k] + deltaZ * uniformDistribution();

	double dx = deltaX / 2;
	double dy = deltaY / 2;
	double dz = deltaZ / 2;

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

	resetElectroMagneticParameters();

	for(int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];

		updateElectroMagneticParameters(particle);
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				velocity[i][j][k] = velocity[i][j][k]/((protonConcentration[i][j][k]*massProton + electronConcentration[i][j][k]*massProton)*volume(i, j, k));
			}
		}
	}

	updateBoundariesParameters();

	//double fullChargeDensity = evaluateFullChargeDensity();
}

void Simulation::updateElectricFluxX(Particle* particle) {
	int xcount = truncate(particle->coordinates.x/deltaX);
	int ycount = truncate((particle->coordinates.y/deltaY) + 0.5);
	int zcount = truncate((particle->coordinates.z/deltaZ) + 0.5);

	double flux = particle->charge*particle->weight*particle->velocityX(speed_of_light_normalized);

	for(int i = xcount - 1; i <= xcount + 1; ++i) {
		for(int j = ycount - 1; j <= ycount + 1; ++j) {
			for(int k = zcount - 1; k <= zcount + 1; ++k) {
				double correlation = correlationWithExBin(i, j, k, particle);
				addElectricFluxX(i, j, k, correlation*flux);
			}
		}
	}
}

void Simulation::updateElectricFluxY(Particle* particle) {
	int xcount = truncate((particle->coordinates.x/deltaX) + 0.5);
	int ycount = truncate((particle->coordinates.y/deltaY));
	int zcount = truncate((particle->coordinates.z/deltaZ) + 0.5);

	double flux = particle->charge*particle->weight*particle->velocityY(speed_of_light_normalized);

	for(int i = xcount - 1; i <= xcount + 1; ++i) {
		for(int j = ycount - 1; j <= ycount + 1; ++j) {
			for(int k = zcount - 1; k <= zcount + 1; ++k) {
				double correlation = correlationWithEyBin(i, j, k, particle);
				addElectricFluxY(i, j, k, correlation*flux);
			}
		}
	}
}

void Simulation::updateElectricFluxZ(Particle* particle) {
	int xcount = truncate((particle->coordinates.x/deltaX) + 0.5);
	int ycount = truncate((particle->coordinates.y/deltaY) + 0.5);
	int zcount = truncate((particle->coordinates.z/deltaZ));

	double flux = particle->charge*particle->weight*particle->velocityZ(speed_of_light_normalized);

	for(int i = xcount - 1; i <= xcount + 1; ++i) {
		for(int j = ycount - 1; j <= ycount + 1; ++j) {
			for(int k = zcount - 1; k <= zcount + 1; ++k) {
				double correlation = correlationWithEzBin(i, j, k, particle);
				addElectricFluxZ(i, j, k, correlation*flux);
			}
		}
	}
}

void Simulation::updateChargeDensity(Particle* particle) {
	int xcount = truncate((particle->coordinates.x/deltaX) + 0.5);
	int ycount = truncate((particle->coordinates.y/deltaY) + 0.5);
	int zcount = truncate((particle->coordinates.z/deltaZ) + 0.5);

	if(xcount < 0) {
		printf("xcount < 0\n");
		exit(0);
	}
	if(xcount > xnumber) {
		printf("xcount > xnumber\n");
		exit(0);
	}

	if(ycount < 0) {
		printf("ycount < 0\n");
		exit(0);
	}
	if(ycount > ynumber) {
		printf("ycount > ynumber\n");
		exit(0);
	}

	if(zcount < 0) {
		printf("zcount < 0\n");
		exit(0);
	}
	if(zcount > znumber) {
		printf("zcount > znumber\n");
		exit(0);
	}

	double charge = particle->charge*particle->weight;

	double fullCorrelation = 0;
	for(int i = xcount - 1; i <= xcount + 1; ++i) {
		for(int j = ycount - 1; j <= ycount + 1; ++j) {
			for(int k = zcount - 1; k <= zcount + 1; ++k) {
				double correlation = correlationWithShiftGridBin(i, j, k, particle);
				fullCorrelation += correlation;
				addChargeDensity(i, j, k, correlation*charge);
				addConcentration(i, j, k, correlation*particle->weight, particle->type);
				addVelocity(i, j, k, particle->momentum*correlation*particle->weight);
			}
		}
	}
}

void Simulation::updateBoundariesParameters() {
	if(boundaryConditionType == PERIODIC){
		for(int j = 0; j < ynumber+1; ++j) {
			for(int k = 0; k < znumber+1; ++k) {
				if(j < ynumber){
					electricFluxY[xnumber][j][k] += electricFluxY[0][j][k];
					//electricFluxY[xnumber][j][k] /=2;
					electricFluxY[0][j][k] = electricFluxY[xnumber][j][k];
				}

				if(k < znumber){
					electricFluxZ[xnumber][j][k] += electricFluxZ[0][j][k];
					//electricFluxZ[xnumber][j][k] /= 2;
					electricFluxZ[0][j][k] = electricFluxZ[znumber][j][k];
				}

				chargeDensity[xnumber][j][k] += chargeDensity[0][j][k];
				//chargeDensity[xnumber][j][k] /= 2;
				chargeDensity[0][j][k] = chargeDensity[xnumber][j][k];

				electronConcentration[xnumber][j][k] += electronConcentration[0][j][k];
				//electronConcentration[xnumber][j][k] /= 2;
				electronConcentration[0][j][k] = electronConcentration[xnumber][j][k];

				protonConcentration[xnumber][j][k] += protonConcentration[0][j][k];
				//protonConcentration[xnumber][j][k] /= 2;
				protonConcentration[0][j][k] = protonConcentration[xnumber][j][k];
			}
		}
	}
	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			if(i < xnumber){
				electricFluxX[i][j][0] += electricFluxX[i][j][znumber];
				electricFluxX[i][j][znumber] = electricFluxX[i][j][0];
			}

			if(j < ynumber){
				electricFluxY[i][j][0] += electricFluxY[i][j][znumber];
				electricFluxY[i][j][znumber] = electricFluxY[i][j][znumber];
			}

			chargeDensity[i][j][0] += chargeDensity[i][j][znumber];
			chargeDensity[i][j][znumber] = chargeDensity[i][j][0];

			electronConcentration[i][j][0] += electronConcentration[i][j][znumber];
			electronConcentration[i][j][znumber] = electronConcentration[i][j][0];

			protonConcentration[i][j][0] += protonConcentration[i][j][znumber];
			protonConcentration[i][j][znumber] = protonConcentration[i][j][0];

		}

		for(int k = 0; k < znumber + 1; ++k) {
			if(i < xnumber) {
				electricFluxX[i][0][k] += electricFluxX[i][ynumber][k];
				electricFluxX[i][ynumber][k] = electricFluxX[i][0][k];
			}

			if(k < znumber){
				electricFluxZ[i][0][k] += electricFluxZ[i][ynumber][k];
				electricFluxZ[i][ynumber][k] = electricFluxZ[i][0][k];
			}

			chargeDensity[i][0][k] += chargeDensity[i][ynumber][k];
			chargeDensity[i][ynumber][k] = chargeDensity[i][0][k];

			electronConcentration[i][0][k] += electronConcentration[i][ynumber][k];
			electronConcentration[i][ynumber][k] = electronConcentration[i][0][k];

			protonConcentration[i][0][k] += protonConcentration[i][ynumber][k];
			protonConcentration[i][ynumber][k] = protonConcentration[i][0][k];
		}
	}
}

void Simulation::updateElectroMagneticParameters(Particle* particle) {
	Particle* tempParticle = new Particle(*particle);

	tempParticle->coordinates = particle->oldCoordinates;
	tempParticle->weight /= 2;
	particle->weight /= 2;

	updateElectricFluxX(tempParticle);
	updateElectricFluxY(tempParticle);
	updateElectricFluxZ(tempParticle);

	updateElectricFluxX(particle);
	updateElectricFluxY(particle);
	updateElectricFluxZ(particle);

	delete tempParticle;

	particle->weight *= 2;
	updateChargeDensity(particle);
}

void Simulation::addElectricFluxX(int i, int j, int k, double flux) {
	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		if(i < 0) {
			i = 0;
			flux = -flux;
		}
		if(i >= xnumber) {
			return;
		}
	}  else if(boundaryConditionType == PERIODIC){
		if(i < 0) {
			i = xnumber - 1;
		}
		if(i >= xnumber) {
			i = 0;
		}
	}

	if(j < 0) {
		j = ynumber - 1;
	}
	if(j > ynumber) {
		j = 1;
	}

	if(k < 0) {
		k = znumber - 1;
	}
	if(k > znumber) {
		k = 1;
	}

	electricFluxX[i][j][k] += flux/volume(i, j, k);

	alertNaNOrInfinity(electricFluxX[i][j][k], "electric flux x = NaN\n");
}

void Simulation::addElectricFluxY(int i, int j, int k, double flux) {
	double factor = 1;
	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		if(i < 0) {
			i = 0;
			//	todo!
			return;
		}

		if(i == 0) {
			factor *= 2;
		}

		if(i == xnumber) {
			factor *= 2;
		}

		if(i > xnumber) {
			return;
		}
	} else if(boundaryConditionType == PERIODIC){
		if(i < 0) {
			i = xnumber - 1;
		}
		if(i > xnumber) {
			i = 1;
		}
	}

	if(j < 0) {
		j = ynumber - 1;
	}
	if(j >= ynumber) {
		j = 0;
	}

	if(k < 0) {
		k = znumber - 1;
	}
	if(k > znumber) {
		k = 1;
	}

	electricFluxY[i][j][k] += factor*flux/volume(i, j, k);

	alertNaNOrInfinity(electricFluxY[i][j][k], "electric flux y = NaN\n");
}

void Simulation::addElectricFluxZ(int i, int j, int k, double flux) {
	double factor = 1;
	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		if(i < 0) {
			i = 0;
			//	todo!
			return;
		}

		if(i == 0) {
			factor *= 2;
		}

		if(i == xnumber) {
			factor *= 2;
		}

		if(i > xnumber) {
			return;
		}
	} else if(boundaryConditionType == PERIODIC){
		if(i < 0) {
			i = xnumber - 1;
		}
		if(i > xnumber) {
			i = 1;
		}
	}

	if(j < 0) {
		j = ynumber - 1;
	}
	if(j > ynumber) {
		j = 1;
	}

	if(k < 0) {
		k = znumber - 1;
	}
	if(k >= znumber) {
		k = 0;
	}

	electricFluxZ[i][j][k] += factor*flux/volume(i, j, k);

	alertNaNOrInfinity(electricFluxZ[i][j][k], "electric flux z = NaN\n");
}

void Simulation::addChargeDensity(int i, int j, int k, double charge) {
	double factor = 1;
	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		if(i < 0) {
			i = 0;
			//	todo!
			return;
		}

		if(i == 0) {
			factor *= 2;
		}

		if(i == xnumber) {
			factor *= 2;
		}

		if(i > xnumber) {
			return;
		}
	} else if(boundaryConditionType == PERIODIC){
		if(i < 0) {
			i = xnumber - 1;
		}
		if(i > xnumber) {
			i = 1;
		}
	}

	if(j < 0) {
		j = ynumber - 1;
	}
	if(j > ynumber) {
		j = 1;
	}

	if(k < 0) {
		k = znumber - 1;
	}
	if(k > znumber) {
		k = 1;
	}

	chargeDensity[i][j][k] += factor*charge/volume(i, j, k);

	alertNaNOrInfinity(chargeDensity[i][j][k], "charge density = NaN\n");	
}

void Simulation::addVelocity(int i, int j, int k, Vector3d momentum) {
	double factor = 1;
	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		if(i < 0) {
			i = 0;
			//	todo!
			return;
		}

		if(i == 0) {
			factor *= 2;
		}

		if(i == xnumber) {
			factor *= 2;
		}

		if(i > xnumber) {
			return;
		}
	} else if(boundaryConditionType == PERIODIC){
		if(i < 0) {
			i = xnumber - 1;
		}
		if(i > xnumber) {
			i = 1;
		}
	}

	if(j < 0) {
		j = ynumber - 1;
	}
	if(j > ynumber) {
		j = 1;
	}

	if(k < 0) {
		k = znumber - 1;
	}
	if(k > znumber) {
		k = 1;
	}

	velocity[i][j][k] += momentum;

	alertNaNOrInfinity(velocity[i][j][k].x, "velocity.x = NaN\n");

	alertNaNOrInfinity(velocity[i][j][k].y, "velocity.y = NaN\n");

	alertNaNOrInfinity(velocity[i][j][k].z, "velocity.z = NaN\n");
}

void Simulation::addConcentration(int i, int j, int k, double weight, ParticleTypes particle_type) {
	double factor = 1;
	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		if(i < 0) {
			i = 0;
			//	todo!
			return;
		}

		if(i == 0) {
			factor *= 2;
		}

		if(i == xnumber) {
			factor *= 2;
		}

		if(i > xnumber) {
			return;
		}
	} else if(boundaryConditionType == PERIODIC){
		if(i < 0) {
			i = xnumber - 1;
		}
		if(i > xnumber) {
			i = 1;
		}
	}

	if(j < 0) {
		j = ynumber - 1;
	}
	if(j > ynumber) {
		j = 1;
	}

	if(k < 0) {
		k = znumber - 1;
	}
	if(k > znumber) {
		k = 1;
	}

	switch (particle_type){
		case ELECTRON :
			electronConcentration[i][j][k] += factor*weight/volume(i, j, k);
			break;
		case PROTON :
			protonConcentration[i][j][k] += factor*weight/volume(i, j, k);
			break;
	}
}
double Simulation::evaluateFullChargeDensity() {
	double fullDensity = 0;
	double fullElectronConcentration = 0;
	double fullProtonConcentration = 0;

	for(int j = 0; j < ynumber; ++j) {
		for(int k = 0; k < znumber; ++k) {
			fullDensity += chargeDensity[0][j][k]*volume(0, j, k)/2;
			fullDensity += chargeDensity[xnumber][j][k]*volume(xnumber, j, k)/2;

			fullElectronConcentration += electronConcentration[0][j][k]*volume(0, j, k)/2;
			fullElectronConcentration += electronConcentration[xnumber][j][k]*volume(xnumber, j, k)/2;

			fullProtonConcentration += protonConcentration[0][j][k]*volume(0, j, k)/2;
			fullProtonConcentration += protonConcentration[xnumber][j][k]*volume(xnumber, j, k)/2;
		}
	}

	for(int i = 1; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				fullDensity += chargeDensity[i][j][k]*volume(i, j, k);

				fullElectronConcentration += electronConcentration[i][j][k]*volume(i, j, k);

				fullProtonConcentration += protonConcentration[i][j][k]*volume(i, j, k);
			}
		}
	}

	fullDensity /= (xsize*ysize*zsize);
	fullElectronConcentration /= (xsize*ysize*zsize);
	fullProtonConcentration /= (xsize*ysize*zsize);

	double fullChargeDensity2 = 0;
	double fullElectronConcentration2 = 0;
	double fullProtonconcentration2 = 0;
	for(int pcount = 0; pcount < particles.size(); ++pcount) {
		fullChargeDensity2 += particles[pcount]->charge*particles[pcount]->weight;
		if(particles[pcount]->type == ELECTRON){
			fullElectronConcentration2 += particles[pcount]->weight;
		}
		if(particles[pcount]->type == PROTON) {
			fullProtonconcentration2 += particles[pcount]->weight;
		}
	}
	fullChargeDensity2 /= (xsize*ysize*zsize);
	fullElectronConcentration2 /= (xsize*ysize*zsize);
	fullProtonconcentration2 /= (xsize*ysize*zsize);

	return fullDensity;
}

Particle* Simulation::electronMinMomentum() {
	Particle* minElectron = getFirstElectron();
	double minMomentum = minElectron->momentum.norm();
	for(int pcount = 0; pcount < particles.size(); ++pcount) {
		Particle* particle = particles[pcount];
		if(particle->type == ELECTRON){
			double momentum = particle->momentum.norm();
			if(momentum < minMomentum) {
				minMomentum = momentum;
				minElectron = particle;
			}
		}
	}

	return minElectron;
}

void Simulation::resetElectroMagneticParameters() {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				electricFluxX[i][j][k] = 0;
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				electricFluxY[i][j][k] = 0;
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber; ++k) {
				electricFluxZ[i][j][k] = 0;
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				chargeDensity[i][j][k] = 0;
				electronConcentration[i][j][k] = 0;
				protonConcentration[i][j][k] = 0;
				velocity[i][j][k] = Vector3d(0, 0, 0);
			}
		}
	}
}

void Simulation::updateEnergy() {
	particleEnergy = 0;
	electricFieldEnergy = 0;
	magneticFieldEnergy = 0;
	energy = 0;

	momentum = Vector3d(0, 0, 0);
	for(int i = 0; i < xnumber+1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				if(i < xnumber) {
					electricFieldEnergy += EfieldX[i][j][k]*EfieldX[i][j][k]*volume(i, j, k);
				}
				double factor = 1;
				if(i == 0 || i == xnumber) {
					factor = 0.5;
				}
				electricFieldEnergy += EfieldY[i][j][k]*EfieldY[i][j][k]*volume(i, j, k)*factor;
				electricFieldEnergy += EfieldZ[i][j][k]*EfieldZ[i][j][k]*volume(i, j, k)*factor;
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				double factor = 1;
				if(i == 0 || i == xnumber) {
					factor = 0.5;
				}
				magneticFieldEnergy += BfieldX[i][j][k]*BfieldX[i][j][k]*volume(i, j, k)*factor;

				if(i < xnumber) {
					magneticFieldEnergy += BfieldY[i][j][k]*BfieldY[i][j][k]*volume(i, j, k);
					magneticFieldEnergy += BfieldZ[i][j][k]*BfieldZ[i][j][k]*volume(i, j, k);
				}
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				//not in one point
				Vector3d E = Vector3d(EfieldX[i][j][k], EfieldY[i][j][k], EfieldZ[i][j][k]);
				Vector3d B = Vector3d(BfieldX[i][j][k], BfieldY[i][j][k], BfieldZ[i][j][k]);

				momentum += E.vectorMult(B)/(4*pi*speed_of_light_normalized);
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