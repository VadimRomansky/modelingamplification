#include "simulation.h"
#include "util.h"
#include "constants.h"

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

	generalizedMinimalResidualMethod(divergenceCleanUpMatrix, divergenceCleanUpRightPart, divergenceCleaningEfield);

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k){
				Efield[i][j][k] += divergenceCleaningEfield[i][j][k];
			}
		}
	}

	updateBoundariesOldField();
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

	divergenceCleanUpRightPart[i][j][k][0] = cleanUpRightPart(i,j, k);
	divergenceCleanUpRightPart[i][j][k][1] = 0;
	divergenceCleanUpRightPart[i][j][k][2] = 0;

	MatrixElement element = MatrixElement(-0.25/deltaX, i, j, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, k, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, j, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-0.25/deltaX, i, nextJ, nextK, 0);
	divergenceCleanUpMatrix[i][j][k][0].push_back(element);

	if(i < xnumber-1){
		element = MatrixElement(0.25/deltaX, i+1, j, k, 0);
		divergenceCleanUpMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(0.25/deltaX, i+1, nextJ, k, 0);
		divergenceCleanUpMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(0.25/deltaX, i+1, j, nextK, 0);
		divergenceCleanUpMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(0.25/deltaX, i+1, nextJ, nextK, 0);
		divergenceCleanUpMatrix[i][j][k][0].push_back(element);
	}


	divergenceCleanUpMatrix[i][j][k][1].push_back(MatrixElement(1, i, j, k, 1));
	divergenceCleanUpMatrix[i][j][k][2].push_back(MatrixElement(1, i, j, k, 2));
}

double Simulation::cleanUpRightPart(int i, int j, int k) {
	double ErightX = (Efield[i+1][j][k].x + Efield[i+1][j+1][k].x + Efield[i+1][j][k+1].x + Efield[i+1][j+1][k+1].x)/4;
	double EleftX = (Efield[i][j][k].x + Efield[i][j][k].x + Efield[i][j][k+1].x + Efield[i][j+1][k+1].x)/4;

	double ErightY = (Efield[i][j+1][k].y + Efield[i+1][j+1][k].y + Efield[i][j+1][k+1].y + Efield[i+1][j+1][k+1].y)/4;
	double EleftY = (Efield[i][j][k].y + Efield[i+1][j][k].y + Efield[i][j][k+1].y + Efield[i+1][j][k+1].y)/4;

	double ErightZ = (Efield[i][j][k+1].z + Efield[i+1][j][k+1].z + Efield[i][j+1][k+1].z + Efield[i+1][j+1][k+1].z)/4;
	double EleftZ = (Efield[i][j][k].z + Efield[i+1][j][k].z + Efield[i][j+1][k].z + Efield[i+1][j+1][k].z)/4;

	double div = ((ErightX - EleftX)/deltaX) + ((ErightY - EleftY)/deltaY) + ((ErightZ - EleftZ)/deltaZ);

	return 4*pi*chargeDensity[i][j][k] - div;
}