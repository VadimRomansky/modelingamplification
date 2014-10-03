#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "specialmath.h"

void Simulation::cleanupDivergence() {
	printf("cleaning up divergence\n");

	for (int i = 0; i <= xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				divergenceCleanUpMatrix[i][j][k][0].clear();
				divergenceCleanUpRightPart[i][j][k][0] = 0;
				createDivergenceCleanupEquation(i, j, k);
			}
		}
	}

	/*if(debugMode) {
		checkEquationMatrix(divergenceCleanUpMatrix);
	}*/

	generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential, xnumber+1, ynumber, znumber, 1);

	updateFieldByPotential();

	updateBoundariesOldField();
}

void Simulation::updateFieldByPotential() {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				Efield[i][j][k] -= evaluateGradPotential(i, j, k);
			}
		}
	}
}

void Simulation::createDivergenceCleanupEquation(int i, int j, int k) {

	if(i == 0) {
		divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, k, 0));
		divergenceCleanUpRightPart[i][j][k][0] = 0;
		return;
	}

	int prevJ = j - 1;
	if(prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int prevK = k - 1;
	if(prevK < 0) {
		prevK = znumber - 1;
	}
	int nextJ = j + 1;
	if(nextJ >= ynumber) {
		nextJ = 0;
	}
	int nextK = k + 1;
	if(nextK >= znumber) {
		nextK = 0;
	}

	double dx2 = deltaX*deltaX;
	double dy2 = deltaY*deltaY;
	double dz2 = deltaZ*deltaZ;

	MatrixElement element;
	if(i == xnumber){
		element = MatrixElement((-1/dx2) - (2/dy2) - (2/dz2), i, j, k, 0);
		divergenceCleanUpMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(1/dx2, i-1, j, k, 0);
		divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	} else if(i == 1){
		element = MatrixElement((-4/dx2) -(2/dy2) - (2/dz2), i, j, k, 0);
		divergenceCleanUpMatrix[i][j][k][0].push_back(element);
		element = MatrixElement((8/(3*dx2)), i-1, j, k, 0);
		divergenceCleanUpMatrix[i][j][k][0].push_back(element);
		element = MatrixElement((4/(3*dx2)), i+1, j, k, 0);
		divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	} else {
		element = MatrixElement((-2/dx2) - (2/dy2) - (2/dz2), i, j, k, 0);
		divergenceCleanUpMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(1/dx2, i-1, j, k, 0);
		divergenceCleanUpMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(1/dx2, i+1, j, k, 0);
		divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	}
	element = MatrixElement(1/dy2, i, prevJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(1/dy2, i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(1/dz2,i, j, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(1/dz2,i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	divergenceCleanUpRightPart[i][j][k][0] = cleanUpRightPart(i, j, k);
}

double Simulation::cleanUpRightPart(int i, int j, int k) {
	double div = evaluateDivE(i-1, j, k);

	return -4*pi*chargeDensity[i-1][j][k] + div;
}