#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "specialmath.h"

void Simulation::cleanupDivergence() {
	printf("cleaning up divergence\n");

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for(int l = 0; l < 3; ++l){
					divergenceCleanUpMatrix[i][j][k][l].clear();
					divergenceCleanUpRightPart[i][j][k][l] = 0;
				}
				createDivergenceCleanupEquation(i, j, k);
			}
		}
	}

	if(debugMode) {
		checkEquationMatrix(divergenceCleanUpMatrix);
	}

	generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningPotential, xnumber, ynumber, znumber, 1);

	updateFieldByPotential();

	updateBoundariesOldField();
}

void Simulation::updateFieldByPotential() {
	for(int i = 1; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				Efield[i][j][k] -= evaluateGradPotential(i, j, k);
			}
		}
	}
}

void Simulation::createDivergenceCleanupEquation(int i, int j, int k) {

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

}

double Simulation::cleanUpRightPart(int i, int j, int k) {
	double div = evaluateDivE(i, j, k);

	return -4*pi*chargeDensity[i][j][k] + div;
}