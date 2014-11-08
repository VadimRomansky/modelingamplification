#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "specialmath.h"


void Simulation::cleanupDivergence() {
	printf("cleaning up divergence\n");

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				divergenceCleanUpMatrix[i][j][k][0].clear();
				divergenceCleanUpMatrix[i][j][k][1].clear();
				divergenceCleanUpMatrix[i][j][k][2].clear();
				divergenceCleanUpRightPart[i][j][k][0] = 0;
				divergenceCleanUpRightPart[i][j][k][1] = 0;
				divergenceCleanUpRightPart[i][j][k][2] = 0;
				if(i == 0 && boundaryConditionType == SUPERCONDUCTERLEFT){
					createDivergenceCleanupLeftEquation(j, k);
				} else if((i == xnumber - 1) && (boundaryConditionType == SUPERCONDUCTERLEFT)) {
					createDivergenceCleanupRightEquation(j, k);
				} else {
					createDivergenceCleanupInternalEquation(i, j, k);
				}
			}
		}
	}

	double fullDensity = 0;
	double fullDiv = 0;
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fullDensity += chargeDensity[i][j][k];
				fullDiv += evaluateDivE(i, j, k);
			}
		}
	}

	int matrixDimension = xnumber*ynumber*znumber;
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				chargeDensity[i][j][k] -= (fullDensity/matrixDimension);
			}
		}
	}

	if(debugMode) {
		checkEquationMatrix(divergenceCleanUpMatrix);
	}

	generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningField, xnumber, ynumber, znumber, 3);

	double**** leftPart = multiplySpecialMatrixVector(divergenceCleanUpMatrix, divergenceCleaningField, xnumber, ynumber, znumber, 3);
	double div = evaluateDivCleaningE(1, 1, 1);
	updateFieldByCleaning();

	updateBoundariesOldField();
}

void Simulation::updateFieldByCleaning() {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				Efield[i][j][k].x += divergenceCleaningField[i][j][k][0];
				Efield[i][j][k].y += divergenceCleaningField[i][j][k][1];
				Efield[i][j][k].z += divergenceCleaningField[i][j][k][2];
			}
		}
	}
}

void Simulation::createDivergenceCleanupInternalEquation(int i, int j, int k) {

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

	int nextI = i+1;
	if(nextI >= xnumber) {
		nextI = 0;
	}

	//div for x
	MatrixElement element = MatrixElement( -0.25/deltaX, i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaX, nextI, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaX, nextI, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaX, nextI, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaX, nextI, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-0.25/deltaY, i, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaY, i, j, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaY, nextI, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaY, nextI, j, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaY, i, nextJ, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaY, i, nextJ, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaY, nextI, nextJ, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaY, nextI, nextJ, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-0.25/deltaZ, i, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaZ, i, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaZ, nextI, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaZ, nextI, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaZ, i, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaZ, i, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaZ, nextI, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaZ, nextI, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	divergenceCleanUpRightPart[i][j][k][0] = cleanUpRightPart(i, j, k);

	//rot z for y

	element = MatrixElement(-0.25/deltaX, i, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, j, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, k, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);

	element = MatrixElement(0.25/deltaX, nextI, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(0.25/deltaX, nextI, j, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(0.25/deltaX, nextI, nextJ, k, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(0.25/deltaX, nextI, nextJ, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);

	element = MatrixElement(-0.25/deltaY, i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-0.25/deltaY, i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-0.25/deltaY, nextI, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-0.25/deltaY, nextI, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);

	element = MatrixElement(0.25/deltaY, i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(0.25/deltaY, i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(0.25/deltaY, nextI, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(0.25/deltaY, nextI, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);


	divergenceCleanUpRightPart[i][j][k][1] = 0;

	//rot y for z;

	element = MatrixElement(-0.25/deltaZ, i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-0.25/deltaZ, i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-0.25/deltaZ, nextI, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-0.25/deltaZ, nextI, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);

	element = MatrixElement(0.25/deltaZ, i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaZ, i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaZ, nextI, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaZ, nextI, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);

	element = MatrixElement(0.25/deltaX, i, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaX, i, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaX, i, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaX, i, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);

	element = MatrixElement(-0.25/deltaX, nextI, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-0.25/deltaX, nextI, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-0.25/deltaX, nextI, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-0.25/deltaX, nextI, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);

	divergenceCleanUpRightPart[i][j][k][2] = 0;
}

void Simulation::createDivergenceCleanupLeftEquation(int j, int k) {
	int i = 0;

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

	int nextI = i+1;

	//div for x

	MatrixElement element = MatrixElement( -0.25/deltaX, i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaX, nextI, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaX, nextI, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaX, nextI, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaX, nextI, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-0.25/deltaY, i, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaY, i, j, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaY, nextI, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaY, nextI, j, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaY, i, nextJ, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaY, i, nextJ, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaY, nextI, nextJ, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaY, nextI, nextJ, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-0.25/deltaZ, i, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaZ, i, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaZ, nextI, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaZ, nextI, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaZ, i, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaZ, i, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaZ, nextI, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaZ, nextI, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	divergenceCleanUpRightPart[i][j][k][0] = cleanUpRightPart(i, j, k);

	//for y

	element = MatrixElement(1, i, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);

	divergenceCleanUpRightPart[i][j][k][1] = 0;

	//for z;

	element = MatrixElement(1, i, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);

	divergenceCleanUpRightPart[i][j][k][2] = 0;
}

void Simulation::createDivergenceCleanupRightEquation(int j, int k) {
	int i = xnumber - 1;

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

	//div for x
	MatrixElement element = MatrixElement( -0.25/deltaX, i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-0.25/deltaY, i, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaY, i, j, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaY, i, nextJ, k, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaY, i, nextJ, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-0.25/deltaZ, i, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaZ, i, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(0.25/deltaZ, i, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(0.25/deltaZ, i, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	divergenceCleanUpRightPart[i][j][k][0] = cleanUpRightPart(i, j, k);

	//rot x for y

	element = MatrixElement(-0.25/deltaY, i, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-0.25/deltaY, i, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);

	element = MatrixElement(0.25/deltaY, i, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(0.25/deltaY, i, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);

	element = MatrixElement(-0.25/deltaZ, i, j, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-0.25/deltaZ, i, nextJ, nextK, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);

	element = MatrixElement(0.25/deltaZ, i, j, k, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(0.25/deltaZ, i, nextJ, k, 1);
	divergenceCleanUpMatrix[i][j][k][1].push_back(element);

	divergenceCleanUpRightPart[i][j][k][1] = 0;

	//rot y for z;

	element = MatrixElement(-0.25/deltaZ, i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-0.25/deltaZ, i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);

	element = MatrixElement(0.25/deltaZ, i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaZ, i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);

	element = MatrixElement(0.25/deltaX, i, j, k, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaX, i, nextJ, k, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaX, i, j, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(0.25/deltaX, i, nextJ, nextK, 2);
	divergenceCleanUpMatrix[i][j][k][2].push_back(element);

	divergenceCleanUpRightPart[i][j][k][2] = 0;
}

double Simulation::cleanUpRightPart(int i, int j, int k) {
	double div = evaluateDivE(i, j, k);

	return 4*pi*chargeDensity[i][j][k] - div;
}