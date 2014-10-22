#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "specialmath.h"

void Simulation::cleanupDivergence() {
	printf("cleaning up divergence\n");

	for (int i = 0; i <= xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				divergenceCleanUpMatrix[i][j][k].clear();
				divergenceCleanUpRightPart[i][j][k] = 0;
				if(boundaryConditionType == SUPERCONDUCTERLEFT){
					if(i == 0) {
						divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(1.0, i, j, k));
						divergenceCleanUpRightPart[i][j][k] = 0;				
					} else if(i == xnumber) {
						createDivergenceCleanupRightEquation(j, k);
					} else {
						createDivergenceCleanupInternalEquation(i, j, k);
					}
				} else if(boundaryConditionType == PERIODIC) {
					if(i == xnumber) {
						divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(1.0, i, j, k));
						divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(-1.0, 0, j, k));
						divergenceCleanUpRightPart[i][j][k] = 0;
					} else if(i == 0 && j == 0 && k == 0) {
						divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(1.0, i, j, k));
						divergenceCleanUpRightPart[i][j][k] = 0;
					} else if(i == xnumber - 1){
						createDivergenceCleanupRightPeriodicEquation(j, k);
					} else {
						createDivergenceCleanupInternalEquation(i,j, k);
					}
				}
			}
		}
	}

	/*if(debugMode) {
		checkEquationMatrix(divergenceCleanUpMatrix);
	}*/

	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential, xnumber+1, ynumber, znumber);
	} if(boundaryConditionType == PERIODIC) {
		generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential, xnumber, ynumber, znumber);
	}

	if(boundaryConditionType == PERIODIC) {
		for(int j = 0; j < ynumber+1; ++j) {
			for(int k = 0; k < znumber+1; ++k) {
				divergenceCleaningPotential[xnumber][j][k] = divergenceCleaningPotential[0][j][k];
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			divergenceCleaningPotential[i][j][znumber] = divergenceCleaningPotential[i][j][0];
		}

		for(int k = 0; k < znumber + 1; ++k) {
			divergenceCleaningPotential[i][ynumber][k] = divergenceCleaningPotential[i][0][k];
		}
	}

	updateFieldByPotential();
}

void Simulation::updateFieldByPotential() {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				EfieldX[i][j][k] -= (divergenceCleaningPotential[i+1][j][k] - divergenceCleaningPotential[i][j][k])/deltaX;
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber + 1; ++k){
				EfieldY[i][j][k] -= (divergenceCleaningPotential[i][j+1][k] - divergenceCleaningPotential[i][j][k])/deltaY;
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber; ++k) {
				EfieldZ[i][j][k] -= (divergenceCleaningPotential[i][j][k+1] - divergenceCleaningPotential[i][j][k])/deltaZ;
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			if(i < xnumber){
				EfieldX[i][j][znumber] = EfieldX[i][j][0];
			}
			if(j < ynumber){
				EfieldY[i][j][znumber] = EfieldY[i][j][0];
			}
		}

		for(int k = 0; k < znumber + 1; ++k) {
			if(i < xnumber) {
				EfieldX[i][ynumber][k] = EfieldX[i][0][k];
			}
			if(k < znumber){
				EfieldZ[i][ynumber][k] = EfieldZ[i][0][k];
			}
		}
	}
}

void Simulation::createDivergenceCleanupInternalEquation(int i, int j, int k) {

	int prevI = i-1;
	if(prevI < 0) {
		prevI = xnumber - 1;
	}
	int nextI = i+1;
	if(nextI > xnumber) {
		nextI = 1;
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

	double element;
	element = -2/(dx2) - 2/(dy2) - 2/(dz2);
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, j, k));

	element = 1/(dx2);
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, nextI, j, k));
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, prevI, j, k));

	element = 1/(dy2);
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, nextJ, k));
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, prevJ, k));

	element = 1/(dz2);
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, j, nextK));
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, j, prevK));
	

	divergenceCleanUpRightPart[i][j][k] = cleanUpRightPart(i, j, k);
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

	double element;
	element = -1/(dx2) - 2/(dy2) - 2/(dz2);
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, j, k));

	element = 1/(dx2);
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i-1, j, k));

	element = 1/(dy2);
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, nextJ, k));
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, prevJ, k));

	element = 1/(dz2);
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, j, nextK));
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, j, prevK));
	

	divergenceCleanUpRightPart[i][j][k] = cleanUpRightPart(i, j, k);
	//divergenceCleanUpRightPart[i][j][k] -= E0.x/deltaX;
}

void Simulation::createDivergenceCleanupRightPeriodicEquation(int j, int k) {
	int i = xnumber - 1;
	int prevI = i - 1;
	int nextI = 0;

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

	double element;
	element = -2/(dx2) - 2/(dy2) - 2/(dz2);
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, j, k));

	element = 1/(dx2);
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, nextI, j, k));
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, prevI, j, k));

	element = 1/(dy2);
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, nextJ, k));
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, prevJ, k));

	element = 1/(dz2);
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, j, nextK));
	divergenceCleanUpMatrix[i][j][k].push_back(MatrixElement(element, i, j, prevK));
	

	divergenceCleanUpRightPart[i][j][k] = cleanUpRightPart(i, j, k);
}

double Simulation::cleanUpRightPart(int i, int j, int k) {
	double div = evaluateDivE(i, j, k);

	return -4*pi*chargeDensity[i][j][k] + div;
}