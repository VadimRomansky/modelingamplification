#include "stdlib.h"
#include "stdio.h"

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "output.h"
#include "constants.h"
#include "matrix3d.h"
#include "random.h"

Simulation::Simulation(){
	debugMode = true;

	currentIteration = 0;
	time = 0;
	particlesNumber = 0;

	theta = 0.5;

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
	E0 = Vector3d(0, 0, 0);


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
	for(int pcount = 0; pcount < particles.size(); ++pcount){
		Particle* particle = particles[pcount];
		delete particle;
	}
	particles.clear();

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				particlesInBbin[i][j][k].clear();
			}
			delete[] Bfield[i][j];
			delete[] newBfield[i][j];
			delete[] tempBfield[i][j];
			delete[] electricDensity[i][j];
			delete[] pressureTensor[i][j];

		}
		delete[] Bfield[i];
		delete[] newBfield[i];
		delete[] tempBfield[i];
		delete[] electricDensity[i];
		delete[] pressureTensor[i];
	}

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber + 1; ++k){
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

	for(int i = 0; i <= xnumber; ++i){
		for(int j = 0; j <= ynumber; ++j){
			for(int k = 0; k <= znumber; ++k){
				for(int l = 0; l < 3; ++l){
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

	delete[] Efield;
	delete[] newEfield;
	delete[] tempEfield;
	delete[] Bfield;
	delete[] newBfield;
	delete[] tempBfield;

	delete[] electricDensity;
	delete[] electricFlux;
	delete[] dielectricTensor;
	delete[] pressureTensor;

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

	for(int i = 0; i < xnumber+1; ++i){
		for(int j = 0; j < ynumber+1; ++j){
			for(int k = 0; k < znumber+1; ++k){
				Efield[i][j][k] = E0;
				newEfield[i][j][k] = Efield[i][j][k];
				tempEfield[i][j][k] = Efield[i][j][k];

			}
		}
	}

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				Bfield[i][j][k] = B0;
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

	Efield = new Vector3d**[xnumber + 1];
	newEfield = new Vector3d**[xnumber + 1];
	tempEfield = new Vector3d**[xnumber + 1];
	Bfield = new Vector3d**[xnumber];
	newBfield = new Vector3d**[xnumber];
	tempBfield = new Vector3d**[xnumber];

	electricFlux = new Vector3d**[xnumber + 1];
	electricDensity = new double**[xnumber];
	dielectricTensor = new Matrix3d**[xnumber + 1];
	pressureTensor = new Matrix3d**[xnumber];


	for(int i = 0; i < xnumber; ++i){
		Bfield[i] = new Vector3d*[ynumber];
		newBfield[i] = new Vector3d*[ynumber];
		tempBfield[i] = new Vector3d*[ynumber];
		electricDensity[i] = new double*[ynumber];
		pressureTensor[i] = new Matrix3d*[ynumber];

		for(int j = 0; j < ynumber; ++j){
			Bfield[i][j] = new Vector3d[znumber];
			newBfield[i][j] = new Vector3d[znumber];
			tempBfield[i][j] = new Vector3d[znumber];
			electricDensity[i][j] = new double[znumber];
			pressureTensor[i][j] = new Matrix3d[znumber];
		}
	}

	for(int i = 0; i < xnumber+1; ++i){
		Efield[i] = new Vector3d*[ynumber+1];
		newEfield[i] = new Vector3d*[ynumber+1];
		tempEfield[i] = new Vector3d*[ynumber+1];
		electricFlux[i] = new Vector3d*[ynumber+1];
		dielectricTensor[i] = new Matrix3d*[ynumber+1];

		for(int j = 0; j < ynumber+1; ++j){
			Efield[i][j] = new Vector3d[znumber+1];
			newEfield[i][j] = new Vector3d[znumber+1];
			tempEfield[i][j] = new Vector3d[znumber+1];
			electricFlux[i][j] = new Vector3d[znumber+1];
			dielectricTensor[i][j] = new Matrix3d[znumber+1];
		}
	}

	maxwellEquationMatrix = new std::vector<MatrixElement>***[xnumber + 1];
	maxwellEquationRightPart = new double***[xnumber + 1];
	for(int i = 0; i <= xnumber; ++i){
		maxwellEquationMatrix[i] = new std::vector<MatrixElement>**[ynumber + 1];
		maxwellEquationRightPart[i] = new double**[ynumber + 1];
		for(int j = 0; j <= ynumber; ++j){
			maxwellEquationMatrix[i][j] = new std::vector<MatrixElement>*[znumber + 1];
			maxwellEquationRightPart[i][j] = new double*[znumber + 1];
			for(int k = 0; k <= znumber; ++k){
				maxwellEquationMatrix[i][j][k] = new std::vector<MatrixElement>[3];
				maxwellEquationRightPart[i][j][k] = new double[3];
			}
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

		evaluateParticlesRotationTensor();
		//evaluateFields();
		moveParticles();
		//updateFields();

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
	return correlationTempEfield(*particle);
}

Vector3d Simulation::correlationBfield(Particle* particle){
	return correlationBfield(*particle);
}

Vector3d Simulation::correlationEfield(Particle* particle){
	return correlationEfield(*particle);
}

Vector3d Simulation::correlationTempEfield(Particle& particle){
	checkParticleInBox(particle);

	int xcount = trunc(particle.coordinates.x/deltaX + 0.5);
	int ycount = trunc(particle.coordinates.y/deltaY + 0.5);
	int zcount = trunc(particle.coordinates.z/deltaZ + 0.5);

	if(xcount < 0){
		printf("xcount < 0\n");
		//exit(0);
	}

	if(xcount > xnumber){
		printf("xcount > xnumber\n");
		//exit(0);
	}

	if(ycount < 0){
		printf("ycount < 0\n");
		//exit(0);
	}

	if(ycount > ynumber){
		printf("ycount > ynumber\n");
		//exit(0);
	}

	if(zcount < 0){
		printf("zcount < 0\n");
		//exit(0);
	}

	if(zcount > znumber){
		printf("zcount > znumber\n");
		//exit(0);
	}

	Vector3d result = Vector3d(0,0,0);

	for(int i = xcount-1; i <= xcount + 1; ++i){
		for(int j = ycount-1; j <= ycount + 1; ++j){
			for(int k = zcount-1; k <= zcount + 1; ++k){
				result = result + correlationFieldWithTempEbin(particle, i, j, k);
			}
		}
	}

	return result;
}

Vector3d Simulation::correlationBfield(Particle& particle){
	checkParticleInBox(particle);

	int xcount = trunc(particle.coordinates.x/deltaX);
	int ycount = trunc(particle.coordinates.y/deltaY);
	int zcount = trunc(particle.coordinates.z/deltaZ);

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

	for(int i = xcount-1; i <= xcount + 1; ++i){
		for(int j = ycount-1; j <= ycount + 1; ++j){
			for(int k = zcount-1; k <= zcount + 1; ++k){
				result = result + correlationFieldWithBbin(particle, i, j, k);
			}
		}
	}

	return result;
}

Vector3d Simulation::correlationEfield(Particle& particle){
	checkParticleInBox(particle);

	int xcount = trunc(particle.coordinates.x/deltaX + 0.5);
	int ycount = trunc(particle.coordinates.y/deltaY + 0.5);
	int zcount = trunc(particle.coordinates.z/deltaZ + 0.5);

	if(xcount < 0){
		printf("xcount < 0\n");
		//exit(0);
	}

	if(xcount > xnumber){
		printf("xcount > xnumber\n");
		//exit(0);
	}

	if(ycount < 0){
		printf("ycount < 0\n");
		//exit(0);
	}

	if(ycount > ynumber){
		printf("ycount > ynumber\n");
		//exit(0);
	}

	if(zcount < 0){
		printf("zcount < 0\n");
		//exit(0);
	}

	if(zcount > znumber){
		printf("zcount > znumber\n");
		//exit(0);
	}

	Vector3d result = Vector3d(0,0,0);

	for(int i = xcount-1; i <= xcount + 1; ++i){
		for(int j = ycount-1; j <= ycount + 1; ++j){
			for(int k = zcount-1; k <= zcount + 1; ++k){
				result = result + correlationFieldWithEbin(particle, i, j, k);
			}
		}
	}

	return result;
}

	

Vector3d Simulation::correlationFieldWithBbin(Particle& particle, int i, int j, int k){

	Vector3d _field;
	
	if(i < 0){
		//_field = Vector3d(0,0,0);
		//note: not zero because reflection
		if(j < 0){
			if(k < 0){
				_field = Bfield[0][ynumber - 1][znumber - 1];
			} else if(k >= znumber){
				_field = Bfield[0][ynumber - 1][0];
			} else {
				_field = Bfield[0][ynumber - 1][k];
			}
		} else if(j >= ynumber){
			if(k < 0){
				_field = Bfield[0][0][znumber - 1];
			} else if(k >= znumber){
				_field = Bfield[0][0][0];
			} else {
				_field = Bfield[0][0][k];
			}
		} else {
			if(k < 0){
				_field = Bfield[0][j][znumber - 1];
			} else if(k >= znumber){
				_field = Bfield[0][j][0];
			} else {
				_field = Bfield[0][j][k];
			}
		}
	} else if( i >= xnumber){
		_field = B0;
	} else {
		if(j < 0){
			if(k < 0){
				_field = Bfield[i][ynumber - 1][znumber - 1];
			} else if(k >= znumber){
				_field = Bfield[i][ynumber - 1][0];
			} else {
				_field = Bfield[i][ynumber - 1][k];
			}
		} else if(j >= ynumber){
			if(k < 0){
				_field = Bfield[i][0][znumber - 1];
			} else if(k >= znumber){
				_field = Bfield[i][0][0];
			} else {
				_field = Bfield[i][0][k];
			}
		} else {
			if(k < 0){
				_field = Bfield[i][j][znumber - 1];
			} else if(k >= znumber){
				_field = Bfield[i][j][0];
			} else {
				_field = Bfield[i][j][k];
			}
		}
	}

	double correlation = correlationWithBbin(particle, i, j, k);

	return _field*correlation;
}

Vector3d Simulation::correlationFieldWithEbin(Particle& particle, int i, int j, int k){

	Vector3d _field = getEfield(i, j, k);

	double correlation = correlationWithEbin(particle, i, j, k);

	return _field*correlation;
}

Vector3d Simulation::correlationFieldWithTempEbin(Particle& particle, int i, int j, int k){

	Vector3d _field = getTempEfield(i, j, k);

	double correlation = correlationWithEbin(particle, i, j, k);

	return _field*correlation;
}

double Simulation::correlationWithBbin(Particle& particle, int i, int j, int k){
	if(! particleCrossBbin(particle, i, j, k))
		return 0.0;

	double x = particle.coordinates.x;
	double y = particle.coordinates.y;
	double z = particle.coordinates.z;

	double leftx;
	double rightx;
	double lefty;
	double righty;
	double leftz;
	double rightz;
		
	if(i < 0){
		leftx = particle.coordinates.x - 2*deltaX;
		rightx = xgrid[0];
	} else if(i >= xnumber){
		leftx = xgrid[xnumber];
		rightx = particle.coordinates.x + 2*deltaX;
	} else {
		leftx = xgrid[i];
		rightx = xgrid[i+1];
		
	}


	if(j < 0){
		lefty = ygrid[0] - deltaY;
		righty = ygrid[0];
	} else if(j >= ynumber){
		lefty = ygrid[ynumber];
		righty = ygrid[ynumber] + deltaY;
	} else if(j == 0){
		if(particle.coordinates.y - particle.dy > ygrid[1]){
			lefty = ygrid[ynumber];
			righty = ygrid[ynumber] + deltaY;
		} else {
			lefty = ygrid[0];
			righty = ygrid[1];
		}
	} else if(j == ynumber - 1){
		if(particle.coordinates.y + particle.dy < ygrid[ynumber-1]){
			lefty = ygrid[0] - deltaY;
			righty = ygrid[0];
		} else {
			lefty = ygrid[ynumber-1];
			righty = ygrid[ynumber];
		}
	} else {
		lefty = ygrid[j];
		righty = ygrid[j+1];
	}

	if(k < 0){
		leftz = zgrid[0] - deltaZ;
		rightz = zgrid[0];
	} else if(k >= znumber){
		leftz = zgrid[znumber];
		rightz = zgrid[znumber] + deltaZ;
	} else if(k == 0){
		if(particle.coordinates.z - particle.dz > zgrid[1]){
			leftz = zgrid[ynumber];
			rightz = zgrid[ynumber] + deltaZ;
		} else {
			leftz = zgrid[0];
			rightz = zgrid[1];
		}
	} else if(k == znumber-1){
		if(particle.coordinates.z + particle.dz < zgrid[znumber - 1]){
			leftz = zgrid[0] - deltaZ;
			rightz = zgrid[0];
		} else {
			leftz = zgrid[znumber-1];
			rightz = zgrid[znumber];
		}
	} else {
		leftz = zgrid[k];
		rightz = zgrid[k+1];
	}


	double correlation = 1;

	double correlationx = correlationBspline(x, particle.dx, leftx, rightx);
	double correlationy = correlationBspline(y, particle.dy, lefty, righty);
	double correlationz = correlationBspline(z, particle.dz, leftz, rightz);

	correlation = correlationx*correlationy*correlationz;

	return correlation;
}

double Simulation::correlationWithEbin(Particle& particle, int i, int j, int k){
	double x = particle.coordinates.x;
	double y = particle.coordinates.y;
	double z = particle.coordinates.z;

	double leftx;
	double rightx;
	double lefty;
	double righty;
	double leftz;
	double rightz;

	if(i < 0){
		leftx = particle.coordinates.x - 2*deltaX;
		rightx = xgrid[0] - deltaX/2;
	} else if(i > xnumber){
		leftx = xgrid[xnumber] + deltaX/2;
		rightx = particle.coordinates.x + 2*deltaX;
	/*} else if(i == 0){
		////note: needs because in the middle of 0 Ebin is wall
		leftx = 0;
		rightx = deltaX/2;*/
	} else {
		leftx = xgrid[i] - deltaX/2;
		rightx = xgrid[i] + deltaX/2;
	} 


	if(j < 0){
		lefty = particle.coordinates.y - 2*deltaY;
		righty = ygrid[0];
	} else if(j > ynumber){
		lefty = ygrid[ynumber] + deltaX/2;
		righty = particle.coordinates.y + 2*deltaY;
	} else {
		lefty = ygrid[j] - deltaX/2;
		righty = ygrid[j] + deltaX/2;
	} 

	if(k < 0){
		leftz = particle.coordinates.z - 2*deltaZ;
		rightz = zgrid[0];
	} else if(k > znumber){
		leftz = zgrid[znumber] + deltaZ/2;
		rightz = particle.coordinates.z + 2*deltaZ;
	} else {
		leftz = zgrid[k] - deltaZ/2;
		rightz = zgrid[k] + deltaZ/2;
	}


	double correlation = 1;

	double correlationx = correlationBspline(x, particle.dx, leftx, rightx);
	double correlationy = correlationBspline(y, particle.dy, lefty, righty);
	double correlationz = correlationBspline(z, particle.dz, leftz, rightz);

	correlation = correlationx*correlationy*correlationz;

	return correlation;
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

Matrix3d Simulation::evaluateAlphaRotationTensor(double beta, Vector3d velocity, Vector3d EField, Vector3d BField){
	Matrix3d result;

	beta = beta/speed_of_light;

	double gamma_factor = 1/sqrt(1 - (velocity.x*velocity.x + velocity.y*velocity.y + velocity.z*velocity.z)/speed_of_light_sqr);
	double G = (beta*(EField.scalarMult(velocity))/speed_of_light_sqr + gamma_factor);
	beta = beta/G;
	double denominator = G*(1 + beta*beta*BField.scalarMult(BField));

	double B[3] = {BField.x, BField.y, BField.z};
	
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
	deltaT = 0.01*delta/speed_of_light;
	deltaT = min2(deltaT, 0.05*massElectron*speed_of_light/(electron_charge*B0.getNorm()));
}

void Simulation::createParticles(){
	printf("creating particles\n");
	//for(int i = 0; i < xnumber; ++i){
		//for(int j = 0; j < ynumber; ++j){
			//for(int k = 0; k < znumber; ++k){
	for(int i = 0; i < 1; ++i){
		for(int j = 0; j < 1; ++j){
			for(int k = 0; k < 1; ++k){
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
		energy = maxwellJuttnerDistribution(temperature, mass);
	}

	double p = sqrt(energy*energy - sqr(mass*speed_of_light_sqr))/speed_of_light;
	
	double pz = p*(2*uniformDistribution() - 1);
	double phi = 2*pi*uniformDistribution();
	double pnormal = sqrt(p*p - pz*pz);
	double px = pnormal*cos(phi);
	double py = pnormal*sin(phi);
	pz = fabs(pz);
	py = fabs(py);
	px = fabs(px);

	Particle* particle = new Particle(mass, charge, weight, type, x, y, z, px, py, pz, dx, dy, dz);

	return particle;
}

double Simulation::volume(int i, int j, int k){
	return deltaX*deltaY*deltaZ;
}

void Simulation::collectParticlesIntoBins(){
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				particlesInBbin[i][j][k].clear();
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i){
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				particlesInEbin[i][j][k].clear();
			}
		}
	}

	for(int pcount = 0; pcount < particles.size(); ++pcount){
		Particle* particle = particles[pcount];
		checkParticleInBox(*particle);

		int xcount = trunc(particle->coordinates.x/deltaX);
		int ycount = trunc(particle->coordinates.y/deltaY);
		int zcount = trunc(particle->coordinates.z/deltaZ);

		for(int i = xcount - 1; i <= xcount + 1; ++i){
			for(int j = ycount - 1; j <= ycount + 1; ++j){
				for(int k = zcount - 1; k <= zcount + 1; ++k){
					if( particleCrossBbin(*particle, i, j, k)){
						particlesInBbin[i][j][k].push_back(particle);
					}
				}
			}
		}

		xcount = trunc(particle->coordinates.x/deltaX + 0.5);
		ycount = trunc(particle->coordinates.y/deltaY + 0.5);
		zcount = trunc(particle->coordinates.z/deltaZ + 0.5);

		for(int i = xcount - 1; i <= xcount + 1; ++i){
			for(int j = ycount - 1; j <= ycount + 1; ++j){
				for(int k = zcount - 1; k <= zcount + 1; ++k){
					if( particleCrossBbin(*particle, i, j, k)){
						particlesInEbin[i][j][k].push_back(particle);
					}
				}
			}
		}
	}

}

bool Simulation::particleCrossBbin(Particle& particle, int i, int j, int k){

	if((xgrid[i] > particle.coordinates.x + particle.dx) || (xgrid[i+1] < particle.coordinates.x - particle.dx))
		return false;

	if( j == 0){
		if((ygrid[j+1] < particle.coordinates.y - particle.dy) && (ygrid[ynumber] > particle.coordinates.y + particle.dy))
			return false;
	} else if( j == ynumber - 1){
		if((ygrid[j] > particle.coordinates.y + particle.dy) && (ygrid[0] < particle.coordinates.y - particle.dy))
			return false;
	} else {
		if((ygrid[j] > particle.coordinates.y + particle.dy) || (ygrid[j+1] < particle.coordinates.y - particle.dy))
			return false;
	}

	if( k == 0){
		if((zgrid[k+1] < particle.coordinates.z - particle.dz) && (zgrid[znumber] > particle.coordinates.z + particle.dz))
			return false;
	} else if( k == znumber - 1){
		if((zgrid[k] > particle.coordinates.z + particle.dz) && (zgrid[0] < particle.coordinates.z - particle.dz))
			return false;
	} else {
		if((zgrid[k] > particle.coordinates.z + particle.dz) || (zgrid[k+1] < particle.coordinates.z - particle.dz))
			return false;
	}

	return true;
}

bool Simulation::particleCrossEbin(Particle& particle, int i, int j, int k){
	if((xgrid[i] - deltaX/2 > particle.coordinates.x + particle.dx) || (xgrid[i+1] - deltaX/2 < particle.coordinates.x - particle.dx))
		return false;
	if((ygrid[j] - deltaY/2 > particle.coordinates.y + particle.dy) || (ygrid[j+1] - deltaY/2 < particle.coordinates.y - particle.dy))
		return false;
	if((zgrid[k] - deltaZ/2 > particle.coordinates.z + particle.dz) || (zgrid[k+1] - deltaZ/2 < particle.coordinates.z - particle.dz))
		return false;
	return true;
}

void Simulation::checkParticleInBox(Particle& particle){
	if(particle.coordinates.x < 0){
		printf("particle.x < 0\n");
		exit(0);
	}
	if(particle.coordinates.x > xgrid[xnumber]){
		printf("particle.x > xsize\n");
		exit(0);
	}
	if(particle.coordinates.y < 0){
		printf("particle.y < 0\n");
		exit(0);
	}
	if(particle.coordinates.y > ygrid[ynumber]){
		printf("particle.y > ysize\n");
		exit(0);
	}
	if(particle.coordinates.z < 0){
		printf("particle.z < 0\n");
		exit(0);
	}
	if(particle.coordinates.z > zgrid[znumber]){
		printf("particle.z > zsize\n");
		exit(0);
	}
}

void Simulation::updateParameters(){
	for(int i = 0; i < xnumber + 1; ++i){
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				electricFlux[i][j][k] = Vector3d(0,0,0);
				dielectricTensor[i][j][k] = Matrix3d(0,0,0,0,0,0,0,0,0);
				for(int pcount = 0; pcount < particlesInEbin[i][j][k].size(); ++pcount){
					Particle* particle = particlesInEbin[i][j][k][pcount];
					double correlation = correlationWithEbin(*particle, i, j, k)/volume(i, j, k);
					double gamma = particle->gammaFactor();
					Vector3d velocity = particle->velocity();
					Vector3d rotatedVelocity = particle->rotationTensor*velocity*gamma;

					if(i > 0){
						electricFlux[i][j][j] = electricFlux[i][j][k] + rotatedVelocity*particle->weight*particle->charge*correlation;
					}
					dielectricTensor[i][j][k] = dielectricTensor[i][j][k] - particle->rotationTensor*(theta*deltaT*deltaT*2*pi*particle->charge*particle->charge*correlation/particle->mass); 
				}
			}
		}
	}

	//for periodic conditions we must summ sides parameters
	for(int i = 0; i < xnumber + 1; ++i){
		for(int j = 0; j < ynumber + 1; ++j){
			electricFlux[i][j][0] = electricFlux[i][j][0] + electricFlux[i][j][znumber];
			electricFlux[i][j][znumber] = electricFlux[i][j][0];

			dielectricTensor[i][j][0] = dielectricTensor[i][j][0] + dielectricTensor[i][j][znumber];
			dielectricTensor[i][j][znumber] = dielectricTensor[i][j][0];
		}

		for(int k = 0; k < znumber + 1; ++k){
			electricFlux[i][0][k] = electricFlux[i][0][k] + electricFlux[i][ynumber][k];
			electricFlux[i][ynumber][k] = electricFlux[i][0][k];

			dielectricTensor[i][0][k] = dielectricTensor[i][0][k] + dielectricTensor[i][ynumber][k];
			dielectricTensor[i][ynumber][k] = dielectricTensor[i][0][k];
		}
	}

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				electricDensity[i][j][k] = 0;
				pressureTensor[i][j][k] = Matrix3d(0,0,0,0,0,0,0,0,0);
				for(int pcount = 0; pcount < particlesInBbin[i][j][k].size(); ++pcount){
					Particle* particle = particlesInBbin[i][j][k][pcount];
					double correlation = correlationWithBbin(*particle, i, j, k)/volume(i, j, k);
					double gamma = particle->gammaFactor();
					Vector3d velocity = particle->velocity();
					Vector3d rotatedVelocity = particle->rotationTensor*velocity*gamma;

					electricDensity[i][j][k] += particle->weight*particle->charge*correlation;

					pressureTensor[i][j][k] = rotatedVelocity.tensorMult(rotatedVelocity)*particle->weight*particle->charge*correlation; 
				}
			}
		}
	}

	for(int i = 0; i <= xnumber; ++i) {
		for(int j = 0; j <= ynumber; ++ j) {
			for(int k = 0; k <= znumber; ++k) {
				Vector3d divPressureTensor = evaluateDivPressureTensor(i, j, k);
				electricFlux[i][j][k] = electricFlux[i][j][k] - divPressureTensor*deltaT/2;
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				double divJ = evaluateDivFlux(i, j, k);

				electricDensity[i][j][k] -= deltaT*theta*divJ;
			}
		}
	}
}

Vector3d Simulation::getBfield(int i, int j, int k){
	if(i < 0){
		return Vector3d(0.0, 0.0, 0.0);
	} else if(i >= xnumber){
		return B0;
	}

	if(j < 0){
		j = ynumber - 1;
	} else if(j >= ynumber){
		j = 0;
	}

	if(k < 0){
		k = znumber - 1;
	} else if(k >= znumber){
		k = 0;
	}

	return Bfield[i][j][k];
}

Vector3d Simulation::getTempEfield(int i, int j, int k){
	if(i < 0){
		return Vector3d(0.0, 0.0, 0.0);
	} else if(i > xnumber){
		return E0;
	}

	if(j < 0){
		j = ynumber-1;
	} else if(j > ynumber){
		j = 1;
	}

	if(k < 0){
		k = znumber-1;
	} else if(k > znumber){
		k = 1;
	}

	return tempEfield[i][j][k];
}

Vector3d Simulation::getEfield(int i, int j, int k) {
	if(i < 0){
		return Vector3d(0.0, 0.0, 0.0);
	} else if(i > xnumber){
		return E0;
	}

	if(j < 0){
		j = ynumber-1;
	} else if(j > ynumber){
		j = 1;
	}

	if(k < 0){
		k = znumber-1;
	} else if(k > znumber){
		k = 1;
	}

	return Efield[i][j][k];
}

Matrix3d Simulation::getPressureTensor(int i, int j, int k){
	if(i < 0){
		return Matrix3d(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	} else if(i >= xnumber){
		return pressureTensor0;
	}

	if(j < 0){
		j = ynumber - 1;
	} else if(j >= ynumber){
		j = 0;
	}

	if(k < 0){
		k = znumber - 1;
	} else if(k >= znumber){
		k = 0;
	}

	return pressureTensor[i][j][k];
}

double Simulation::getDensity(int i, int j, int k) {
	if(i < 0) return 0;
	if (i >= xnumber) return 0;

	if(j < 0) {
		j = ynumber - 1;
	}
	if(j >= ynumber){
		j = 0;
	}

	if(k < 0){
		k = znumber - 1;
	} 
	if(k >= znumber){
		k = 0;
	}

	return electricDensity[i][j][k];
}