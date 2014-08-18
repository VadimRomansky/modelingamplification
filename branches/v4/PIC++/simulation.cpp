#include "stdlib.h"
#include "stdio.h"

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "output.h"
#include "constants.h"
#include "matrix3d.h"

Simulation::Simulation(){
	currentIteration = 0;
	time = 0;
	particlesNumber = 0;

	//read input!
	xnumber = 10;
	ynumber = 10;
	znumber = 10;

	xsize = 1E9;
	ysize = 1E9;
	zsize = 1E9;

	temperature = 1000;
	density = 1.6E-24;

	maxIteration = 1E7;
	maxTime = 1E7;

	particlesPerBin = 10;


	Kronecker = Matrix3d(1.0,0,0,0,1.0,0,0,0,1.0);

	for(int i = 0; i < 3; ++i){
		for(int j = 0; j < 3; ++j){
			for(int k = 0; k < 3; ++k){
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

Simulation::~Simulation(){
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			delete[] Efield[i][j];
			delete[] newEfield[i][j];
			delete[] tempEfield[i][j];
			delete[] Bfield[i][j];
			delete[] newBfield[i][j];
			delete[] tempBfield[i][j];
		}
		delete[] Efield[i];
		delete[] newEfield[i];
		delete[] tempEfield[i];
		delete[] Bfield[i];
		delete[] newBfield[i];
		delete[] tempBfield[i];
	}

	delete[] Efield;
	delete[] newEfield;
	delete[] tempEfield;
	delete[] Bfield;
	delete[] newBfield;
	delete[] tempBfield;

	delete[] xgrid;
	delete[] ygrid;
	delete[] zgrid;
	delete[] middleXgrid;
	delete[] middleYgrid;
	delete[] middleZgrid;
}

void Simulation::initialize(){

	deltaX = xsize/(xnumber-1);
	deltaY = ysize/(ynumber-1);
	deltaZ = zsize/(znumber-1);

	for(int i = 0; i <= xnumber; ++i){
		xgrid[i] = i*deltaX;
	}

	for(int i = 0; i <= ynumber; ++i){
		ygrid[i] = i*deltaY;
	}

	for(int i = 0; i <= znumber; ++i){
		zgrid[i] = i*deltaZ;
	}

	for(int i = 0; i < xnumber; ++i){
		middleXgrid[i] = (xgrid[i] + xgrid[i+1])/2;
	}

	for(int i = 0; i < ynumber; ++i){
		middleYgrid[i] = (ygrid[i] + ygrid[i+1])/2;
	}

	for(int i = 0; i < znumber; ++i){
		middleZgrid[i] = (zgrid[i] + zgrid[i+1])/2;
	}

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				Efield[i][j][k] = Vector3d();
				newEfield[i][j][k] = Efield[i][j][k];
				tempEfield[i][j][k] = Efield[i][j][k];
				Bfield[i][j][k] = Vector3d(1E-6,0,0);
				newBfield[i][j][k] = Bfield[i][j][k];
				tempBfield[i][j][k] = Bfield[i][j][k];

			}
		}
	}
}

void Simulation::createArrays(){
	xgrid = new double[xnumber + 1];
	ygrid = new double[ynumber + 1];
	zgrid = new double[znumber + 1];

	middleXgrid = new double[xnumber];
	middleYgrid = new double[ynumber];
	middleZgrid = new double[znumber];

	Efield = new Vector3d**[xnumber];
	newEfield = new Vector3d**[xnumber];
	tempEfield = new Vector3d**[xnumber];
	Bfield = new Vector3d**[xnumber];
	newBfield = new Vector3d**[xnumber];
	tempBfield = new Vector3d**[xnumber];
	for(int i = 0; i < xnumber; ++i){
		Efield[i] = new Vector3d*[ynumber];
		newEfield[i] = new Vector3d*[ynumber];
		tempEfield[i] = new Vector3d*[ynumber];
		Bfield[i] = new Vector3d*[ynumber];
		newBfield[i] = new Vector3d*[ynumber];
		tempBfield[i] = new Vector3d*[ynumber];

		for(int j = 0; j < ynumber; ++j){
			Efield[i][j] = new Vector3d[znumber];
			newEfield[i][j] = new Vector3d[znumber];
			tempEfield[i][j] = new Vector3d[znumber];
			Bfield[i][j] = new Vector3d[znumber];
			newBfield[i][j] = new Vector3d[znumber];
			tempBfield[i][j] = new Vector3d[znumber];
		}
	}
}

void Simulation::createFiles(){
	traectoryFile = fopen("./output/traectory.dat","w");
	fclose(traectoryFile);
	distributionFile = fopen("./output/distribution_protons.dat","w");
	fclose(distributionFile);
}

void Simulation::simulate(){
	createArrays();
	initialize();
	createFiles();
	createParticles();

	updateDeltaT();
	//deltaT = 0;

	while(time < maxTime && currentIteration < maxIteration){
		printf("iteration number %d time = %lf\n", currentIteration, time);
		moveParticles();
		evaluateFields();

		time += deltaT;
		currentIteration++;

		if(currentIteration % writeParameter == 0){
			distributionFile = fopen("./output/distribution_protons.dat", "a");
			outputDistribution(distributionFile, particles, ParticleTypes::PROTON);
			fclose(distributionFile);
			traectoryFile = fopen("./output/traectory.dat","a");
			outputTraectory(traectoryFile, particles[0], time);
			fclose(traectoryFile);
		}
	}
}

Vector3d Simulation::correlationTempEfield(Particle* particle){
	return correlationField(particle, tempEfield);
}

Vector3d Simulation::correlationBfield(Particle* particle){
	return correlationField(particle, Bfield);
}

Vector3d Simulation::correlationTempEfield(Particle& particle){
	return correlationField(particle, tempEfield);
}

Vector3d Simulation::correlationBfield(Particle& particle){
	return correlationField(particle, Bfield);
}

Vector3d Simulation::correlationField(Particle* particle, Vector3d*** field){
	return correlationField(*particle, field);
}
	
Vector3d Simulation::correlationField(Particle& particle, Vector3d*** field){
	int xcount = particle.coordinates.x/deltaX;
	int ycount = particle.coordinates.y/deltaY;
	int zcount = particle.coordinates.z/deltaZ;

	if(xcount < 0){
		printf("xcount < 0\n");
		//exit(0);
	}

	if(xcount >= xnumber){
		printf("xcount >= xnumber\n");
		//exit(0);
	}

	if(ycount < 0){
		printf("ycount < 0\n");
		//exit(0);
	}

	if(ycount >= ynumber){
		printf("ycount >= ynumber\n");
		//exit(0);
	}

	if(zcount < 0){
		printf("zcount < 0\n");
		//exit(0);
	}

	if(zcount >= znumber){
		printf("zcount >= znumber\n");
		//exit(0);
	}

	Vector3d result = Vector3d(0,0,0);

	for(int i = max2(0,xcount-1); i <= min2(xnumber-1, xcount + 1); ++i){
		for(int j = max2(0,ycount-1); j <= min2(ynumber-1, ycount + 1); ++j){
			for(int k = max2(0,zcount-1); k <= min2(znumber-1, zcount + 1); ++k){
				result = result + correlationFieldWithBin(particle, field, i, j, k);
			}
		}
	}

	return result;
}

Vector3d Simulation::correlationFieldWithBin(Particle& particle, Vector3d*** field, int i, int j, int k){
	double x = particle.coordinates.x;
	double y = particle.coordinates.y;
	double z = particle.coordinates.z;

	double leftx = xgrid[i];
	double rightx = xgrid[i+1];
	double lefty = ygrid[j];
	double righty = ygrid[j+1];
	double leftz = zgrid[k];
	double rightz = zgrid[k+1];

	double correlation = 1;

	double correlationx = correlationBspline(x, particle.dx, leftx, rightx);
	double correlationy = correlationBspline(y, particle.dy, lefty, righty);
	double correlationz = correlationBspline(z, particle.dz, leftz, rightz);

	correlation = correlationx*correlationy*correlationz;

	return field[i][j][k]*correlation;
}

double Simulation::correlationBspline(const double& x, const double&  dx, const double& leftx, const double& rightx){

	if(rightx < leftx){
		printf("rightx < leftx\n");
		exit(0);
	}
	if(dx > rightx - leftx){
		printf("dx < rightx - leftx\n");
		exit(0);
	}

	double correlation = 0;

	if(x < leftx - dx)
		return 0;
	if(x > rightx + dx)
		return 0;

	if(x < leftx){
		correlation = sqr(x + dx - leftx)/2;
	} else if ( x > rightx){
		correlation = sqr(rightx - (x - dx))/2;
	} else if ( x < leftx + dx){
		correlation = sqr(dx) - sqr(leftx - (x - dx))/2;
	} else if ( x > rightx - dx){
		correlation = sqr(dx) - sqr(x + dx - rightx)/2;
	} else {
		correlation = dx*dx;
	}

	correlation /= dx*dx;

	return correlation;
}

Matrix3d Simulation::evaluateAlphaRotationTensor(double beta, Vector3d BField){
	Matrix3d result;

	beta = beta/speed_of_light;

	double B[3] = {BField.x, BField.y, BField.z};

	double denominator = 1 + beta*beta*(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
	
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			result.matrix[i][j] = Kronecker.matrix[i][j] + beta*beta*B[i]*B[j];
			for(int k = 0; k < 3; ++k){
				for(int l = 0; l < 3; ++l){
					result.matrix[i][j] -= beta*LeviCivita[j][k][l]*Kronecker.matrix[i][l]*B[k];
				}
			}

			result.matrix[i][j] /= denominator;
		}
	}

	return result;
}

void Simulation::updateDeltaT(){
	double delta = min3(deltaX, deltaY, deltaZ);
	deltaT = delta/speed_of_light;
}

void Simulation::createParticles(){
	printf("creating particles\n");
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				for(int l = 0; l < particlesPerBin; ++l){
					ParticleTypes type;
					if(l % 2 == 0){
						type = ParticleTypes::PROTON;
					} else {
						type = ParticleTypes::ELECTRON;
					}
					double weight = density*2/(massProton*particlesPerBin)*volume(i, j, k);
					Particle* particle = createParticle(i, j, k, weight, type);
					particles.push_back(particle);
					particlesNumber++;
					printf("%d\n", particlesNumber);
				}
			}
		}
	}
}

Particle* Simulation::createParticle(int i, int j, int k, double weight, ParticleTypes type){
	double charge = 0;
	double mass = 0;

	switch(type){
		case ParticleTypes::PROTON :
			mass = massProton;
			charge = electron_charge;
			break;
		case ParticleTypes::ELECTRON :
			mass = massElectron;
			charge = -electron_charge;
			break;
	}

	double x = xgrid[i] + deltaX*uniformDistribution();
	double y = ygrid[i] + deltaY*uniformDistribution();
	double z = zgrid[i] + deltaZ*uniformDistribution();

	double dx = deltaX/4;
	double dy = deltaY/4;
	double dz = deltaZ/4;

	double energy = mass*speed_of_light_sqr;

	double theta = kBoltzman*temperature/(mass*speed_of_light_sqr);

	if(theta < 0.01){
		energy = mass*speed_of_light_sqr + maxwellDistribution(temperature);
	} else {
		energy = maxwellJuttnerDistribution(temperature);
	}

	double p = sqrt(energy*energy - sqr(mass*speed_of_light_sqr))/speed_of_light;
	
	double pz = p*(2*uniformDistribution() - 1);
	double phi = 2*pi*uniformDistribution();
	double pnormal = sqrt(p*p - pz*pz);
	double px = pnormal*cos(phi);
	double py = pnormal*sin(phi);

	Particle* particle = new Particle(mass, charge, weight, type, x, y, z, px, py, pz, dx, dy, dz);

	return particle;
}

double Simulation::volume(int i, int j, int k){
	return deltaX*deltaY*deltaZ;
}