#include "stdlib.h"
#include "stdio.h"

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"

void Simulation::evaluateFields(){
	updateParameters();

	evaluateMaxwellEquationMatrix();

	generalizedMinimalResidualMethod();

	evaluateMagneticField();

	for(int i = 0; i < xnumber + 1; ++i){
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				newEfield[i][j][k] = (tempEfield[i][j][k] - tempEfield[i][j][k]*(1 - theta))/theta;
			}
		}
	}
}

void Simulation::evaluateMaxwellEquationMatrix(){
}

void Simulation::evaluateMagneticField(){
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				Vector3d rotE = evaluateRotE(i, j, k);
				newBfield[i][j][k] = Bfield[i][j][k] + rotE*deltaT;
			}
		}
	}
}

double Simulation::evaluateDivFlux(int i, int j, int k){
	if(i < 0){
		printf("i < 0\n");
		exit(0);
	}

	if(i >= xnumber){
		printf("x >= xnumber\n");
		exit(0);
	}

	if(j < 0){
		printf("j < 0\n");
		exit(0);
	}

	if(j >= ynumber){
		printf("j >= ynumber\n");
		exit(0);
	}

	if(k < 0){
		printf("k < 0\n");
		exit(0);
	}

	if(k >= znumber){
		printf("k >= znumber\n");
		exit(0);
	}


	double rfluxX = (electricFlux[i+1][j][k].x + electricFlux[i+1][j+1][k].x + electricFlux[i+1][j][k+1].x + electricFlux[i+1][j+1][k+1].x)/4;
	double lfluxX = (electricFlux[i][j][k].x + electricFlux[i][j+1][k].x + electricFlux[i][j][k+1].x + electricFlux[i][j+1][k+1].x)/4;

	double rfluxY = (electricFlux[i][j+1][k].y + electricFlux[i+1][j+1][k].y + electricFlux[i][j+1][k+1].y + electricFlux[i+1][j+1][k+1].y)/4;
	double lfluxY = (electricFlux[i][j][k].y + electricFlux[i+1][j][k].y + electricFlux[i][j][k+1].y + electricFlux[i+1][j][k+1].y)/4;

	double rfluxZ = (electricFlux[i][j][k+1].z + electricFlux[i+1][j][k+1].z + electricFlux[i][j+1][k+1].z + electricFlux[i+1][j+1][k+1].z)/4;
	double lfluxZ = (electricFlux[i][j][k].z + electricFlux[i+1][j][k].z + electricFlux[i][j+1][k].z + electricFlux[i+1][j+1][k].z)/4;

	return ((rfluxX - lfluxX)/deltaX) + ((rfluxY - lfluxY)/deltaY) + ((rfluxZ - lfluxZ)/deltaZ);
}

Vector3d Simulation::evaluateRotB(int i, int j, int k){
	if(i < 0){
		printf("i < 0\n");
		exit(0);
	}

	if(i >= xnumber){
		printf("x >= xnumber\n");
		exit(0);
	}

	if(j < 0){
		printf("j < 0\n");
		exit(0);
	}

	if(j >= ynumber){
		printf("j >= ynumber\n");
		exit(0);
	}

	if(k < 0){
		printf("k < 0\n");
		exit(0);
	}

	if(k >= znumber){
		printf("k >= znumber\n");
		exit(0);
	}

	Vector3d BrightX;
	Vector3d BleftX;
	Vector3d BrightY;
	Vector3d BleftY;
	Vector3d BrightZ;
	Vector3d BleftZ;

	if(i == 0){
	} else {
		BrightX = (getBfield(i,j,k) + getBfield(i, j-1, k) + getBfield(i, j-1, k-1) + getBfield(i, j-1, k-1))/4;
		BleftX = (getBfield(i-1,j,k) + getBfield(i-1, j-1, k) + getBfield(i-1, j-1, k-1) + getBfield(i-1, j-1, k-1))/4;

		BrightY = (getBfield(i,j,k) + getBfield(i-1, j, k) + getBfield(i, j, k-1) + getBfield(i-1, j, k-1))/4;
		BleftY = (getBfield(i,j-1,k) + getBfield(i-1, j-1, k) + getBfield(i, j-1, k-1) + getBfield(i-1, j-1, k-1))/4;

		BrightZ = (getBfield(i,j,k) + getBfield(i-1, j, k) + getBfield(i, j-1, k) + getBfield(i-1, j-1, k))/4;
		BleftZ = (getBfield(i,j,k-1) + getBfield(i-1, j, k-1) + getBfield(i, j-1, k-1) + getBfield(i-1, j-1, k-1))/4;
	}


	double x = 0;
	double y = 0;
	double z = 0;

	x = ((BrightY.z - BleftY.z)/deltaY) - ((BrightZ.y - BleftZ.y)/deltaZ);
	y = ((BrightZ.x - BleftZ.x)/deltaZ) - ((BrightX.z - BleftX.z)/deltaX);
	z = ((BrightX.y - BleftX.y)/deltaX) - ((BrightY.x - BleftY.x)/deltaY);

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotE(int i, int j, int k){
	if(i < 0){
		printf("i < 0\n");
		exit(0);
	}

	if(i > xnumber){
		printf("x >= xnumber\n");
		exit(0);
	}

	if(j < 0){
		printf("j < 0\n");
		exit(0);
	}

	if(j > ynumber){
		printf("j >= ynumber\n");
		exit(0);
	}

	if(k < 0){
		printf("k < 0\n");
		exit(0);
	}

	if(k > znumber){
		printf("k >= znumber\n");
		exit(0);
	}

	double x = 0;
	double y = 0;
	double z = 0;

	Vector3d ErightX = (tempEfield[i+1][j][k] + tempEfield[i+1][j+1][k] + tempEfield[i+1][j][k+1] + tempEfield[i+1][j+1][k+1])/4;
	Vector3d EleftX = (tempEfield[i][j][k] + tempEfield[i][j+1][k] + tempEfield[i][j][k+1] + tempEfield[i][j+1][k+1])/4;

	Vector3d ErightY = (tempEfield[i][j+1][k] + tempEfield[i+1][j+1][k] + tempEfield[i][j+1][k+1] + tempEfield[i+1][j+1][k+1])/4;
	Vector3d EleftY = (tempEfield[i][j][k] + tempEfield[i+1][j][k] + tempEfield[i][j][k+1] + tempEfield[i+1][j][k+1])/4;

	Vector3d ErightZ = (tempEfield[i][j][k+1] + tempEfield[i+1][j][k+1] + tempEfield[i][j+1][k+1] + tempEfield[i+1][j+1][k+1])/4;
	Vector3d EleftZ = (tempEfield[i][j][k] + tempEfield[i+1][j][k] + tempEfield[i][j+1][k] + tempEfield[i+1][j+1][k])/4;

	x = ((ErightY.z - EleftY.z)/deltaY) - ((ErightZ.y - EleftZ.y)/deltaZ);
	y = ((ErightZ.x - EleftZ.x)/deltaZ) - ((ErightX.z - EleftX.z)/deltaX);
	z = ((ErightX.y - EleftX.y)/deltaX) - ((ErightY.x - EleftY.x)/deltaY);

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateDivPressureTensor(int i, int j, int k){
	Vector3d result = Vector3d(0,0,0);

	Matrix3d tensorDerX = (getPressureTensor(i, j, k) + getPressureTensor(i, j-1, k) + getPressureTensor(i, j, k-1) + getPressureTensor(i, j-1, k-1) - getPressureTensor(i-1, j, k) - getPressureTensor(i-1, j-1, k) - getPressureTensor(i-1, j, k-1) - getPressureTensor(i-1, j-1, k-1))/(4*deltaX);
	Matrix3d tensorDerY = (getPressureTensor(i, j, k) + getPressureTensor(i-1, j, k) + getPressureTensor(i, j, k-1) + getPressureTensor(i-1, j, k-1) - getPressureTensor(i, j-1, k) - getPressureTensor(i-1, j-1, k) - getPressureTensor(i, j-1, k-1) - getPressureTensor(i-1, j-1, k-1))/(4*deltaY);
	Matrix3d tensorDerZ = (getPressureTensor(i, j, k) + getPressureTensor(i-1, j, k) + getPressureTensor(i, j-1, k) + getPressureTensor(i-1, j-1, k) - getPressureTensor(i, j, k-1) - getPressureTensor(i-1, j, k-1) - getPressureTensor(i, j-1, k-1) - getPressureTensor(i-1, j-1, k-1))/(4*deltaZ);

	result.x = tensorDerX.matrix[0][0] + tensorDerY.matrix[1][0] + tensorDerZ.matrix[2][0];
	result.y = tensorDerX.matrix[0][1] + tensorDerY.matrix[1][1] + tensorDerZ.matrix[2][1];
	result.z = tensorDerX.matrix[0][2] + tensorDerY.matrix[1][2] + tensorDerZ.matrix[2][2];

	return result;
}