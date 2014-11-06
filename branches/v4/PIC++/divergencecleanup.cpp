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
				if(i == 0) {
					divergenceCleanUpMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, k, 0));
					divergenceCleanUpRightPart[i][j][k][0] = 0;				
				} else if(i == 1) {
					createDivergenceCleanupLeftEquation(j, k);
				} else if(i == xnumber) {
					createDivergenceCleanupRightEquation(j, k);
				} else {
					createDivergenceCleanupInternalEquation(i, j, k);
				}
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

	double dx2 = deltaX*deltaX;
	double dy2 = deltaY*deltaY;
	double dz2 = deltaZ*deltaZ;

	MatrixElement element;
	
	element = MatrixElement(-(1/(2*dx2)) - (1/(2*dy2)) - (1/(2*dz2)), i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(4*dx2)) - (1/(4*dy2)) - (1/(4*dz2)), i+1, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(4*dx2)) - (1/(4*dy2)) - (1/(4*dz2)), i-1, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-(1/(4*dx2)) + (1/(4*dy2)) - (1/(4*dz2)), i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-(1/(4*dx2)) + (1/(4*dy2)) - (1/(4*dz2)), i, prevJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-(1/(4*dx2)) - (1/(4*dy2)) + (1/(4*dz2)), i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-(1/(4*dx2)) - (1/(4*dy2)) + (1/(4*dz2)), i, j, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(8*dx2)) + (1/(8*dy2)) - (1/(4*dz2)), i+1, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dy2)) - (1/(4*dz2)), i+1, prevJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dy2)) - (1/(4*dz2)), i-1, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dy2)) - (1/(4*dz2)), i-1, prevJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(8*dx2)) + (1/(8*dz2)) - (1/(4*dy2)), i+1, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dz2)) - (1/(4*dy2)), i+1, j, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dz2)) - (1/(4*dy2)), i-1, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dz2)) - (1/(4*dy2)), i-1, j, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(8*dy2)) + (1/(8*dz2)) - (1/(4*dx2)), i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dy2)) + (1/(8*dz2)) - (1/(4*dx2)), i, nextJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dy2)) + (1/(8*dz2)) - (1/(4*dx2)), i, prevJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dy2)) + (1/(8*dz2)) - (1/(4*dx2)), i, prevJ, prevJ, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i+1, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i+1, nextJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i+1, prevJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i+1, prevJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i-1, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i-1, nextJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i-1, prevJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i-1, prevJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);


	divergenceCleanUpRightPart[i][j][k][0] = cleanUpRightPart(i, j, k);
}

void Simulation::createDivergenceCleanupLeftEquation(int j, int k) {
	int i = 1;

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

	element = MatrixElement(-(3/(4*dx2)) - (1/(4*dy2)) - (1/(4*dz2)), i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(4*dx2)) - (1/(4*dy2)) - (1/(4*dz2)), i+1, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(2*dx2)) - (1/(2*dy2)) - (1/(2*dz2)), i-1, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-(3/(8*dx2)) + (1/(8*dy2)) - (1/(8*dz2)), i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-(3/(8*dx2)) + (1/(8*dy2)) - (1/(8*dz2)), i, prevJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-(3/(8*dx2)) - (1/(8*dy2)) + (1/(8*dz2)), i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-(3/(8*dx2)) - (1/(8*dy2)) + (1/(8*dz2)), i, j, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(8*dx2)) + (1/(8*dy2)) - (1/(4*dz2)), i+1, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dy2)) - (1/(4*dz2)), i+1, prevJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(4*dx2)) + (1/(4*dy2)) - (1/(2*dz2)), i-1, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(4*dx2)) + (1/(4*dy2)) - (1/(2*dz2)), i-1, prevJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(8*dx2)) + (1/(8*dz2)) - (1/(4*dy2)), i+1, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dz2)) - (1/(4*dy2)), i+1, j, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(4*dx2)) + (1/(4*dz2)) - (1/(2*dy2)), i-1, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(4*dx2)) + (1/(4*dz2)) - (1/(2*dy2)), i-1, j, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(16*dy2)) + (1/(16*dz2)) - (3/(8*dx2)), i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dy2)) + (1/(16*dz2)) - (3/(8*dx2)), i, nextJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dy2)) + (1/(16*dz2)) - (3/(8*dx2)), i, prevJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dy2)) + (1/(16*dz2)) - (3/(8*dx2)), i, prevJ, prevJ, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i+1, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i+1, nextJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i+1, prevJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i+1, prevJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dy2)) + (1/(8*dz2)), i-1, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dy2)) + (1/(8*dz2)), i-1, nextJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dy2)) + (1/(8*dz2)), i-1, prevJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dy2)) + (1/(8*dz2)), i-1, prevJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	divergenceCleanUpRightPart[1][j][k][0] = cleanUpRightPart(1, j, k);
}

void Simulation::createDivergenceCleanupRightEquation(int j, int k) {
	int i = xnumber;
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
	
	element = MatrixElement(-(1/(dx2)) - (1/(2*dy2)) - (1/(2*dz2)), i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(4*dx2)) - (1/(4*dy2)) - (1/(4*dz2)), i-1, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-(1/(8*dx2)) + (1/(4*dy2)) - (1/(4*dz2)), i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-(1/(8*dx2)) + (1/(4*dy2)) - (1/(4*dz2)), i, prevJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement(-(1/(8*dx2)) - (1/(4*dy2)) + (1/(4*dz2)), i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-(1/(8*dx2)) - (1/(4*dy2)) + (1/(4*dz2)), i, j, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(8*dx2)) + (1/(8*dy2)) - (1/(4*dz2)), i-1, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dy2)) - (1/(4*dz2)), i-1, prevJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(8*dx2)) + (1/(8*dz2)) - (1/(4*dy2)), i-1, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dx2)) + (1/(8*dz2)) - (1/(4*dy2)), i-1, j, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(8*dy2)) + (1/(8*dz2)) - (1/(8*dx2)), i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dy2)) + (1/(8*dz2)) - (1/(8*dx2)), i, nextJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dy2)) + (1/(8*dz2)) - (1/(8*dx2)), i, prevJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(8*dy2)) + (1/(8*dz2)) - (1/(8*dx2)), i, prevJ, prevJ, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i-1, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i-1, nextJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i-1, prevJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement((1/(16*dx2)) + (1/(16*dy2)) + (1/(16*dz2)), i-1, prevJ, prevK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);


	divergenceCleanUpRightPart[i][j][k][0] = cleanUpRightPart(i, j, k);


}

double Simulation::cleanUpRightPart(int i, int j, int k) {
	double div = evaluateDivE(i-1, j, k);

	return -4*pi*chargeDensity[i-1][j][k] + div;
}