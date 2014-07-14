#include "stdlib.h"
#include "stdio.h"

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "constants.h"

Simulation::Simulation(){
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
	createArrays();

	deltaX = 1;
	deltaY = 1;
	deltaZ = 1;

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

		for(int j = 0; j < ynumber; ++i){
			Efield[i][j] = new Vector3d[znumber];
			newEfield[i][j] = new Vector3d[znumber];
			tempEfield[i][j] = new Vector3d[znumber];
			Bfield[i][j] = new Vector3d[znumber];
			newBfield[i][j] = new Vector3d[znumber];
			tempBfield[i][j] = new Vector3d[znumber];

			for(int k = 0; k < znumber; ++k){
				Efield[i][j][k] = Vector3d();
				newEfield[i][j][k] = Efield[i][j][k];
				tempEfield[i][j][k] = Efield[i][j][k];
				Bfield[i][j][k] = Vector3d(1,0,0);
				newBfield[i][j][k] = Bfield[i][j][k];
				tempBfield[i][j][k] = Bfield[i][j][k];

			}
		}
	}
}

Vector3d Simulation::correlationTempEfield(Particle* particle){
	return correlationField(particle, tempEfield);
}

Vector3d Simulation::correlationBfield(Particle* particle){
	return correlationField(particle, Bfield);
}
	
Vector3d Simulation::correlationField(Particle* particle, Vector3d*** field){
	int i = particle->coordinates.x/deltaX;
	int j = particle->coordinates.y/deltaY;
	int k = particle->coordinates.z/deltaZ;

	if(i < 0){
		printf("i < 0\n");
		//exit(0);
	}

	if(i >= xnumber){
		printf("i >= xnumber\n");
		//exit(0);
	}

	if(j < 0){
		printf("j < 0\n");
		//exit(0);
	}

	if(j >= ynumber){
		printf("j >= ynumber\n");
		//exit(0);
	}

	if(k < 0){
		printf("k < 0\n");
		//exit(0);
	}

	if(k >= znumber){
		printf("k >= znumber\n");
		//exit(0);
	}


}

Vector3d Simulation::correlationFieldWithBin(Particle* particle, Vector3d*** field, int i, int j, int k){
	double x = particle->coordinates.x;
	double y = particle->coordinates.y;
	double z = particle->coordinates.z;

	double leftx = xgrid[i];
	double rightx = xgrid[i+1];
	double lefty = ygrid[i];
	double righty = ygrid[i+1];
	double leftz = zgrid[i];
	double rightz = zgrid[i+1];

	double correlation = 1;

	correlation *= correlationBspline(x, particle->dx, leftx, rightx);
	correlation *= correlationBspline(y, particle->dy, lefty, righty);
	correlation *= correlationBspline(z, particle->dz, leftz, rightz);

	return field[i][j][k]*correlation;
}

double Simulation::correlationBspline(const double& x, const double&  dx, const double& leftx, const double& rightx){

	if(rightx < leftx){
		printf("rightx < leftx\n");
		exit(0);
	}
	if(dx < rightx - leftx){
		printf("dx < rightx - leftx\n");
		exit(0);
	}

	double correlation = 0;

	if(x < leftx){
		correlation = sqr(dx - (leftx - x))/2;
	} else if ( x > rightx){
		correlation = sqr(dx - (x - rightx))/2;
	} else if ( x < leftx + dx){
		correlation = (sqr(dx) - sqr(dx - (leftx - x)))/2;
	} else if ( x > rightx - dx){
		correlation = (sqr(dx) - sqr(dx - (x - rightx)))/2;
	} else {
		correlation = dx*dx;
	}

	correlation /= dx*dx;

	return correlation;
}