#include "stdlib.h"
#include "stdio.h"

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"

void Simulation::evaluateFields() {
	updateElectroMagneticParameters();

	evaluateMaxwellEquationMatrix();

	generalizedMinimalResidualMethod(maxwellEquationMatrix, maxwellEquationRightPart, tempEfield);

	updateBoundaries();

	evaluateMagneticField();


	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				newEfield[i][j][k] = (tempEfield[i][j][k] - Efield[i][j][k] * (1 - theta)) / theta;
			}
		}
	}
}

void Simulation::updateFields() {
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				Efield[i][j][k] = newEfield[i][j][k];
			}
		}
	}

	for (int i = 0; i < xnumber ; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber ; ++k) {
				Bfield[i][j][k] = newBfield[i][j][k];
			}
		}
	}
}

void Simulation::evaluateMaxwellEquationMatrix() {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					maxwellEquationMatrix[i][j][k][l].clear();
				}
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				double rightPartVector[3];
				if (i == 0) {
					createPerfectConductaryBoundaryCondition(j, k);
				} else {
					createInternalEquation(i, j, k);
				}
			}
		}
	}

	if (debugMode) {
		checkEquationMatrix(maxwellEquationMatrix);
	}
}

void Simulation::checkEquationMatrix(std::vector<MatrixElement>**** matrix) {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						MatrixElement element = matrix[i][j][k][l][m];
						if (element.i < 0) {
							printf("element i < 0\n");
							exit(0);
						}
						if (element.i >= xnumber) {
							printf("element i >= xnumber");
							exit(0);
						}
						if (element.j < 0) {
							printf("element j < 0\n");
							exit(0);
						}
						if (element.j >= ynumber) {
							printf("element j >= ynumber\n");
							exit(0);
						}
						if (element.k < 0) {
							printf("element k < 0\n");
							exit(0);
						}
						if (element.k >= znumber) {
							printf("eement k >= znumber\n");
							exit(0);
						}
						for (int n = m + 1; n < matrix[i][j][k][l].size(); ++n) {
							MatrixElement tempElement = matrix[i][j][k][l][n];

							if (element.equalsIndex(tempElement)) {
								printf("equals indexes\n");
								printf("current = %d %d %d %d\n", i, j, k, l);
								printf("temp = %d %d %d %d\n", element.i, element.j, element.k, element.l);
								exit(0);
							}
						}
					}
				}
			}
		}
	}
}

void Simulation::createPerfectConductaryBoundaryCondition(int j, int k) {
	int i = 0;
	int nextJ = j + 1;
	int nextK = k + 1;
	if (nextJ >= ynumber) {
		nextJ = 0;
	}
	if (nextK >= znumber) {
		nextK = 0;
	}

	maxwellEquationRightPart[i][j][k][0] = 4 * pi * electricDensity[0][j][k] * deltaX * deltaY * deltaZ;

	double element = -deltaZ * deltaY * (1 + dielectricTensor[0][j][k].matrix[0][0]) / 4;
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, 0, j, k, 0));

	element = -deltaZ * deltaY * (1 + dielectricTensor[0][nextJ][k].matrix[0][0]) / 4;
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, 0, nextJ, k, 0));

	element = -deltaZ * deltaY * (1 + dielectricTensor[0][j][nextK].matrix[0][0]) / 4;
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, 0, j, nextK, 0));

	element = -deltaZ * deltaY * (1 + dielectricTensor[0][nextJ][nextK].matrix[0][0]) / 4;
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, 0, nextJ, nextK, 0));

	element = (deltaZ * deltaY * (1 + dielectricTensor[i + 1][j][k].matrix[0][0]) / 4) - (deltaY * deltaX * (dielectricTensor[i + 1][j][k].matrix[2][0]) / 4) - (deltaX * deltaZ * (dielectricTensor[i + 1][j][k].matrix[1][0]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, j, k, 0));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][j][k].matrix[0][1] / 4) - (deltaY * deltaX * dielectricTensor[i + 1][j][k].matrix[2][1] / 4) - (deltaX * deltaZ * (1 + dielectricTensor[i + 1][j][k].matrix[1][1]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, j, k, 1));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][j][k].matrix[0][2] / 4) - (deltaY * deltaX * (1 + dielectricTensor[i + 1][j][k].matrix[2][2]) / 4) - (deltaX * deltaZ * dielectricTensor[i + 1][j][k].matrix[1][2] / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, j, k, 2));

	element = (deltaZ * deltaY * (1 + dielectricTensor[i + 1][nextJ][k].matrix[0][0]) / 4) - (deltaY * deltaX * (dielectricTensor[i + 1][nextJ][k].matrix[2][0]) / 4) + (deltaX * deltaZ * (dielectricTensor[i + 1][nextJ][k].matrix[1][0]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, nextJ, k, 0));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][nextJ][k].matrix[0][1] / 4) - (deltaY * deltaX * dielectricTensor[i + 1][nextJ][k].matrix[2][1] / 4) + (deltaX * deltaZ * (1 + dielectricTensor[i + 1][nextJ][k].matrix[1][1]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, nextJ, k, 1));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][nextJ][k].matrix[0][2] / 4) - (deltaY * deltaX * (1 + dielectricTensor[i + 1][nextJ][k].matrix[2][2]) / 4) + (deltaX * deltaZ * dielectricTensor[i + 1][nextJ][k].matrix[1][2] / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, nextJ, k, 2));

	element = (deltaZ * deltaY * (1 + dielectricTensor[i + 1][j][nextK].matrix[0][0]) / 4) + (deltaY * deltaX * (dielectricTensor[i + 1][j][nextK].matrix[2][0]) / 4) - (deltaX * deltaZ * (dielectricTensor[i + 1][j][nextK].matrix[1][0]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, j, nextK, 0));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][j][nextK].matrix[0][1] / 4) + (deltaY * deltaX * dielectricTensor[i + 1][j][nextK].matrix[2][1] / 4) - (deltaX * deltaZ * (1 + dielectricTensor[i + 1][j][nextK].matrix[1][1]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, j, nextK, 1));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][j][nextK].matrix[0][2] / 4) + (deltaY * deltaX * (1 + dielectricTensor[i + 1][j][nextK].matrix[2][2]) / 4) - (deltaX * deltaZ * dielectricTensor[i + 1][j][nextK].matrix[1][2] / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, j, nextK, 2));

	element = (deltaZ * deltaY * (1 + dielectricTensor[i + 1][nextJ][nextK].matrix[0][0]) / 4) + (deltaY * deltaX * (dielectricTensor[i + 1][nextJ][nextK].matrix[2][0]) / 4) + (deltaX * deltaZ * (dielectricTensor[i + 1][nextJ][nextK].matrix[1][0]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, nextJ, nextK, 0));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][nextJ][nextK].matrix[0][1] / 4) + (deltaY * deltaX * dielectricTensor[i + 1][nextJ][nextK].matrix[2][1] / 4) + (deltaX * deltaZ * (1 + dielectricTensor[i + 1][nextJ][nextK].matrix[1][1]) / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, nextJ, nextK, 1));

	element = (deltaZ * deltaY * dielectricTensor[i + 1][nextJ][nextK].matrix[0][2] / 4) + (deltaY * deltaX * (1 + dielectricTensor[i + 1][nextJ][nextK].matrix[2][2]) / 4) + (deltaX * deltaZ * dielectricTensor[i + 1][nextJ][nextK].matrix[1][2] / 4);
	maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(element, i + 1, nextJ, nextK, 2));


	/*maxwellEquationMatrix[i][j][k][0].push_back(MatrixElement(1.0, i, j, k, 0));
	maxwellEquationRightPart[i][j][k][0] = E0.x;*/
	//Ey and Ez = 0
	maxwellEquationMatrix[i][j][k][1].push_back(MatrixElement(1.0, i, j, k, 1));
	maxwellEquationRightPart[i][j][k][1] = 0;
	maxwellEquationMatrix[i][j][k][2].push_back(MatrixElement(1.0, i, j, k, 2));
	maxwellEquationRightPart[i][j][k][2] = 0;
}

void Simulation::createInternalEquation(int i, int j, int k) {
	Vector3d rightPart = Efield[i][j][k];

	//rightPart = rightPart + (evaluateRotB(i, j, k) - electricFlux[i][j][k] * 4 * pi / speed_of_light_normalized) * speed_of_light_normalized * theta * deltaT;
	rightPart = rightPart - (evaluateRotB(i, j, k) - electricFlux[i][j][k] * 4 * pi / speed_of_light_normalized) * speed_of_light_normalized * theta * deltaT;
	rightPart = rightPart - evaluateGradDensity(i, j, k) * 4 * pi * sqr(speed_of_light_normalized * theta * deltaT);

	createInternalEquationX(i, j, k, rightPart);
	createInternalEquationY(i, j, k, rightPart);
	createInternalEquationZ(i, j, k, rightPart);

	maxwellEquationRightPart[i][j][k][0] = rightPart.x;
	maxwellEquationRightPart[i][j][k][1] = rightPart.y;
	maxwellEquationRightPart[i][j][k][2] = rightPart.z;
}

void Simulation::createInternalEquationX(int i, int j, int k, Vector3d& rightPart) {
	int prevJ = j - 1;
	int prevK = k - 1;
	int nextJ = j + 1;
	int nextK = k + 1;
	if (prevJ < 0) {
		prevJ = ynumber - 1;
	}
	if (prevK < 0) {
		prevK = znumber - 1;
	}
	if (nextJ >= ynumber) {
		nextJ = 0;
	}
	if (nextK >= znumber) {
		nextK = 0;
	}

	double cthetadt2 = sqr(speed_of_light_normalized * theta * deltaT);

	MatrixElement element = MatrixElement(1 + dielectricTensor[i][j][k].matrix[0][0] + cthetadt2 * ((2 * (1 + 0.25 * dielectricTensor[i][j][k].matrix[0][0]) / (deltaX * deltaX)) + (2 / (deltaY * deltaY)) + (2 / (deltaZ * deltaZ))), i, j, k, 0);

	//E i j k
	maxwellEquationMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(dielectricTensor[i][j][k].matrix[0][1] + cthetadt2 * ((0.5 * dielectricTensor[i][j][k].matrix[0][1] / (deltaX * deltaX))), i, j, k, 1);
	maxwellEquationMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(dielectricTensor[i][j][k].matrix[0][2] + cthetadt2 * ((0.5 * dielectricTensor[i][j][k].matrix[0][2] / (deltaX * deltaX))), i, j, k, 2);
	maxwellEquationMatrix[i][j][k][0].push_back(element);

	//E i+1 j k
	if (i < xnumber - 1) {
		element = MatrixElement(cthetadt2 * (-(1.0 + 0.5 * dielectricTensor[i + 1][j][k].matrix[0][0]) / (deltaX * deltaX)), i + 1, j, k, 0);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(cthetadt2 * (-(0.5 * dielectricTensor[i + 1][j][k].matrix[0][1]) / (deltaX * deltaX)), i + 1, j, k, 1);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(cthetadt2 * (-(0.5 * dielectricTensor[i + 1][j][k].matrix[0][2]) / (deltaX * deltaX)), i + 1, j, k, 2);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
	} else {
		rightPart.x += cthetadt2 * ((1.0 + 0.5 * dielectricTensor[i + 1][j][k].matrix[0][0]) / (deltaX * deltaX)) * E0.x;
		rightPart.x += cthetadt2 * ((0.5 * dielectricTensor[i + 1][j][k].matrix[0][1]) / (deltaX * deltaX))*E0.y;
		rightPart.x += cthetadt2 * ((0.5 * dielectricTensor[i + 1][j][k].matrix[0][2]) / (deltaX * deltaX))*E0.z;
	}

	//E i-1 j k
	element = MatrixElement(cthetadt2 * (-(1.0 + 0.5 * dielectricTensor[i - 1][j][k].matrix[0][0]) / (deltaX * deltaX)), i - 1, j, k, 0);
	maxwellEquationMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(cthetadt2 * (-(0.5 * dielectricTensor[i - 1][j][k].matrix[0][1]) / (deltaX * deltaX)), i - 1, j, k, 1);
	maxwellEquationMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(cthetadt2 * (-(0.5 * dielectricTensor[i - 1][j][k].matrix[0][2]) / (deltaX * deltaX)), i - 1, j, k, 2);
	maxwellEquationMatrix[i][j][k][0].push_back(element);

	//Ex i j+1 k
	element = MatrixElement(-cthetadt2 / (deltaY * deltaY), i, nextJ, k, 0);
	maxwellEquationMatrix[i][j][k][0].push_back(element);

	//Ex i j-1 k
	element = MatrixElement(-cthetadt2 / (deltaY * deltaY), i, prevJ, k, 0);
	maxwellEquationMatrix[i][j][k][0].push_back(element);

	//Ex i j k+1
	element = MatrixElement(-cthetadt2 / (deltaZ * deltaZ), i, j, nextK, 0);
	maxwellEquationMatrix[i][j][k][0].push_back(element);

	//Ex i j k-1
	element = MatrixElement(-cthetadt2 / (deltaZ * deltaZ), i, j, prevK, 0);
	maxwellEquationMatrix[i][j][k][0].push_back(element);

	//E i+1 j+1 k
	if (i < xnumber - 1) {
		element = MatrixElement(-cthetadt2 * (dielectricTensor[i + 1][nextJ][k].matrix[1][0] / (4 * deltaX * deltaY)), i + 1, nextJ, k, 0);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(-cthetadt2 * (dielectricTensor[i + 1][nextJ][k].matrix[1][1] / (4 * deltaX * deltaY)), i + 1, nextJ, k, 1);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(-cthetadt2 * (dielectricTensor[i + 1][nextJ][k].matrix[1][2] / (4 * deltaX * deltaY)), i + 1, nextJ, k, 2);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
	} else {
		rightPart.x += cthetadt2 * (dielectricTensor[i + 1][nextJ][k].matrix[1][0] / (4 * deltaX * deltaY))*E0.x;
		rightPart.x += cthetadt2 * (dielectricTensor[i + 1][nextJ][k].matrix[1][1] / (4 * deltaX * deltaY))*E0.y;
		rightPart.x += cthetadt2 * (dielectricTensor[i + 1][nextJ][k].matrix[1][2] / (4 * deltaX * deltaY))*E0.z;
	}

	//E i-1 j-1 k
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i - 1][prevJ][k].matrix[1][0] / (4 * deltaX * deltaY)), i - 1, prevJ, k, 0);
	maxwellEquationMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i - 1][prevJ][k].matrix[1][1] / (4 * deltaX * deltaY)), i - 1, prevJ, k, 1);
	maxwellEquationMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i - 1][prevJ][k].matrix[1][2] / (4 * deltaX * deltaY)), i - 1, prevJ, k, 2);
	maxwellEquationMatrix[i][j][k][0].push_back(element);

	//E i+1 j-1 k
	if (i < xnumber - 1) {
		element = MatrixElement(cthetadt2 * (dielectricTensor[i + 1][prevJ][k].matrix[1][0] / (4 * deltaX * deltaY)), i + 1, prevJ, k, 0);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(cthetadt2 * (dielectricTensor[i + 1][prevJ][k].matrix[1][1] / (4 * deltaX * deltaY)), i + 1, prevJ, k, 1);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(cthetadt2 * (dielectricTensor[i + 1][prevJ][k].matrix[1][2] / (4 * deltaX * deltaY)), i + 1, prevJ, k, 2);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
	} else {
		rightPart.x -= cthetadt2 * (dielectricTensor[i + 1][prevJ][k].matrix[1][0] / (4 * deltaX * deltaY))*E0.x;
		rightPart.x -= cthetadt2 * (dielectricTensor[i + 1][prevJ][k].matrix[1][1] / (4 * deltaX * deltaY))*E0.y;
		rightPart.x -= cthetadt2 * (dielectricTensor[i + 1][prevJ][k].matrix[1][2] / (4 * deltaX * deltaY))*E0.z;
	}

	//E i-1 j+1 k
	element = MatrixElement(cthetadt2 * (dielectricTensor[i - 1][nextJ][k].matrix[1][0] / (4 * deltaX * deltaY)), i - 1, nextJ, k, 0);
	maxwellEquationMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i - 1][nextJ][k].matrix[1][1] / (4 * deltaX * deltaY)), i - 1, nextJ, k, 1);
	maxwellEquationMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i - 1][nextJ][k].matrix[1][2] / (4 * deltaX * deltaY)), i - 1, nextJ, k, 2);
	maxwellEquationMatrix[i][j][k][0].push_back(element);

	//E i+1 j k+1
	if (i < xnumber < 1) {
		element = MatrixElement(-cthetadt2 * (dielectricTensor[i + 1][j][nextK].matrix[2][0] / (4 * deltaX * deltaZ)), i + 1, j, nextK, 0);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(-cthetadt2 * (dielectricTensor[i + 1][j][nextK].matrix[2][1] / (4 * deltaX * deltaZ)), i + 1, j, nextK, 1);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(-cthetadt2 * (dielectricTensor[i + 1][j][nextK].matrix[2][2] / (4 * deltaX * deltaZ)), i + 1, j, nextK, 2);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
	} else {
		rightPart.x += cthetadt2 * (dielectricTensor[i + 1][j][nextK].matrix[2][0] / (4 * deltaX * deltaZ))*E0.x;
		rightPart.x += cthetadt2 * (dielectricTensor[i + 1][j][nextK].matrix[2][1] / (4 * deltaX * deltaZ))*E0.y;
		rightPart.x += cthetadt2 * (dielectricTensor[i + 1][j][nextK].matrix[2][2] / (4 * deltaX * deltaZ))*E0.z;
	}

	//E i-1 j k-1
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i - 1][j][prevK].matrix[2][0] / (4 * deltaX * deltaZ)), i - 1, j, prevK, 0);
	maxwellEquationMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i - 1][j][prevK].matrix[2][1] / (4 * deltaX * deltaZ)), i - 1, j, prevK, 1);
	maxwellEquationMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i - 1][j][prevK].matrix[2][2] / (4 * deltaX * deltaZ)), i - 1, j, prevK, 2);
	maxwellEquationMatrix[i][j][k][0].push_back(element);

	//E i+1 j k-1
	if (i < xnumber - 1) {
		element = MatrixElement(cthetadt2 * (dielectricTensor[i + 1][j][prevK].matrix[2][0] / (4 * deltaX * deltaZ)), i + 1, j, prevK, 0);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(cthetadt2 * (dielectricTensor[i + 1][j][prevK].matrix[2][1] / (4 * deltaX * deltaZ)), i + 1, j, prevK, 1);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
		element = MatrixElement(cthetadt2 * (dielectricTensor[i + 1][j][prevK].matrix[2][2] / (4 * deltaX * deltaZ)), i + 1, j, prevK, 2);
		maxwellEquationMatrix[i][j][k][0].push_back(element);
	} else {
		rightPart.x -= cthetadt2 * (dielectricTensor[i + 1][j][prevK].matrix[2][0] / (4 * deltaX * deltaZ))*E0.x;
		rightPart.x -= cthetadt2 * (dielectricTensor[i + 1][j][prevK].matrix[2][1] / (4 * deltaX * deltaZ))*E0.y;
		rightPart.x -= cthetadt2 * (dielectricTensor[i + 1][j][prevK].matrix[2][2] / (4 * deltaX * deltaZ))*E0.z;
	}

	//E i-1 j k+1
	element = MatrixElement(cthetadt2 * (dielectricTensor[i - 1][j][nextK].matrix[2][0] / (4 * deltaX * deltaZ)), i - 1, j, nextK, 0);
	maxwellEquationMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i - 1][j][nextK].matrix[2][1] / (4 * deltaX * deltaZ)), i - 1, j, nextK, 1);
	maxwellEquationMatrix[i][j][k][0].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i - 1][j][nextK].matrix[2][2] / (4 * deltaX * deltaZ)), i - 1, j, nextK, 2);
	maxwellEquationMatrix[i][j][k][0].push_back(element);

	double value = 0;
	for(int m = 0; m < maxwellEquationMatrix[i][j][k][0].size(); ++m) {
		value += maxwellEquationMatrix[i][j][k][0][m].value*Efield[maxwellEquationMatrix[i][j][k][0][m].i][maxwellEquationMatrix[i][j][k][0][m].j][maxwellEquationMatrix[i][j][k][0][m].k][maxwellEquationMatrix[i][j][k][0][m].l];
	}
	value = value - rightPart.x;
}

void Simulation::createInternalEquationY(int i, int j, int k, Vector3d& rightPart) {
	int prevJ = j - 1;
	int prevK = k - 1;
	int nextJ = j + 1;
	int nextK = k + 1;
	if (prevJ < 0) {
		prevJ = ynumber - 1;
	}
	if (prevK < 0) {
		prevK = znumber - 1;
	}
	if (nextJ >= ynumber) {
		nextJ = 0;
	}
	if (nextK >= znumber) {
		nextK = 0;
	}
	double cthetadt2 = sqr(speed_of_light_normalized * theta * deltaT);

	MatrixElement element = MatrixElement(1 + dielectricTensor[i][j][k].matrix[1][1] + cthetadt2 * ((2 / (deltaX * deltaX)) + (2 * (1 + 0.25 * dielectricTensor[i][j][k].matrix[1][1]) / (deltaY * deltaY)) + (2 / (deltaZ * deltaZ))), i, j, k, 1);

	//E i j k
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(dielectricTensor[i][j][k].matrix[1][0] + cthetadt2 * ((0.5 * dielectricTensor[i][j][k].matrix[1][0] / (deltaY * deltaY))), i, j, k, 0);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(dielectricTensor[i][j][k].matrix[1][2] + cthetadt2 * ((0.5 * dielectricTensor[i][j][k].matrix[1][2] / (deltaY * deltaY))), i, j, k, 2);
	maxwellEquationMatrix[i][j][k][1].push_back(element);

	//E i j+1 k
	element = MatrixElement(cthetadt2 * (-(0.5 * dielectricTensor[i][nextJ][k].matrix[1][0]) / (deltaY * deltaY)), i, nextJ, k, 0);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(cthetadt2 * (-(1.0 + 0.5 * dielectricTensor[i][nextJ][k].matrix[1][1]) / (deltaY * deltaY)), i, nextJ, k, 1);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(cthetadt2 * (-(0.5 * dielectricTensor[i][nextJ][k].matrix[1][2]) / (deltaY * deltaY)), i, nextJ, k, 2);
	maxwellEquationMatrix[i][j][k][1].push_back(element);

	//E i j-1 k
	element = MatrixElement(cthetadt2 * (-(0.5 * dielectricTensor[i][prevJ][k].matrix[1][0]) / (deltaY * deltaY)), i, prevJ, k, 0);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(cthetadt2 * (-(1.0 + 0.5 * dielectricTensor[i][prevJ][k].matrix[1][1]) / (deltaY * deltaY)), i, prevJ, k, 1);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(cthetadt2 * (-(0.5 * dielectricTensor[i][prevJ][k].matrix[1][2]) / (deltaY * deltaY)), i, prevJ, k, 2);
	maxwellEquationMatrix[i][j][k][1].push_back(element);

	//Ey i+1 j k
	if (i < xnumber - 1) {
		element = MatrixElement(-cthetadt2 / (deltaX * deltaX), i + 1, j, k, 1);
		maxwellEquationMatrix[i][j][k][1].push_back(element);
	} else {
		rightPart.y += (cthetadt2 / (deltaX * deltaX))*E0.y;
	}

	//Ey i-1 j k
	element = MatrixElement(-cthetadt2 / (deltaX * deltaX), i - 1, j, k, 1);
	maxwellEquationMatrix[i][j][k][1].push_back(element);

	//Ey i j k+1
	element = MatrixElement(-cthetadt2 / (deltaZ * deltaZ), i, j, nextK, 1);
	maxwellEquationMatrix[i][j][k][1].push_back(element);

	//Ey i j k-1
	element = MatrixElement(-cthetadt2 / (deltaZ * deltaZ), i, j, prevK, 1);
	maxwellEquationMatrix[i][j][k][1].push_back(element);

	//E i+1 j+1 k
	if (i < xnumber - 1) {
		element = MatrixElement(-cthetadt2 * (dielectricTensor[i + 1][nextJ][k].matrix[0][0] / (4 * deltaX * deltaY)), i + 1, nextJ, k, 0);
		maxwellEquationMatrix[i][j][k][1].push_back(element);
		element = MatrixElement(-cthetadt2 * (dielectricTensor[i + 1][nextJ][k].matrix[0][1] / (4 * deltaX * deltaY)), i + 1, nextJ, k, 1);
		maxwellEquationMatrix[i][j][k][1].push_back(element);
		element = MatrixElement(-cthetadt2 * (dielectricTensor[i + 1][nextJ][k].matrix[0][2] / (4 * deltaX * deltaY)), i + 1, nextJ, k, 2);
		maxwellEquationMatrix[i][j][k][1].push_back(element);
	} else {
		rightPart.y += cthetadt2 * (dielectricTensor[i + 1][nextJ][k].matrix[0][0] / (4 * deltaX * deltaY))*E0.x;
		rightPart.y += cthetadt2 * (dielectricTensor[i + 1][nextJ][k].matrix[0][1] / (4 * deltaX * deltaY))*E0.y;
		rightPart.y += cthetadt2 * (dielectricTensor[i + 1][nextJ][k].matrix[0][2] / (4 * deltaX * deltaY))*E0.z;
	}

	//E i-1 j-1 k
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i - 1][prevJ][k].matrix[0][0] / (4 * deltaX * deltaY)), i - 1, prevJ, k, 0);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i - 1][prevJ][k].matrix[0][1] / (4 * deltaX * deltaY)), i - 1, prevJ, k, 1);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i - 1][prevJ][k].matrix[0][2] / (4 * deltaX * deltaY)), i - 1, prevJ, k, 2);
	maxwellEquationMatrix[i][j][k][1].push_back(element);

	//E i+1 j-1 k
	if (i < xnumber - 1) {
		element = MatrixElement(cthetadt2 * (dielectricTensor[i + 1][prevJ][k].matrix[0][0] / (4 * deltaX * deltaY)), i + 1, prevJ, k, 0);
		maxwellEquationMatrix[i][j][k][1].push_back(element);
		element = MatrixElement(cthetadt2 * (dielectricTensor[i + 1][prevJ][k].matrix[0][1] / (4 * deltaX * deltaY)), i + 1, prevJ, k, 1);
		maxwellEquationMatrix[i][j][k][1].push_back(element);
		element = MatrixElement(cthetadt2 * (dielectricTensor[i + 1][prevJ][k].matrix[0][2] / (4 * deltaX * deltaY)), i + 1, prevJ, k, 2);
		maxwellEquationMatrix[i][j][k][1].push_back(element);
	} else {
		rightPart.y -= cthetadt2 * (dielectricTensor[i + 1][prevJ][k].matrix[0][0] / (4 * deltaX * deltaY))*E0.x;
		rightPart.y -= cthetadt2 * (dielectricTensor[i + 1][prevJ][k].matrix[0][1] / (4 * deltaX * deltaY))*E0.y;
		rightPart.y -= cthetadt2 * (dielectricTensor[i + 1][prevJ][k].matrix[0][2] / (4 * deltaX * deltaY))*E0.z;
	}

	//E i-1 j+1 k
	element = MatrixElement(cthetadt2 * (dielectricTensor[i - 1][nextJ][k].matrix[0][0] / (4 * deltaX * deltaY)), i - 1, nextJ, k, 0);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i - 1][nextJ][k].matrix[0][1] / (4 * deltaX * deltaY)), i - 1, nextJ, k, 1);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i - 1][nextJ][k].matrix[0][2] / (4 * deltaX * deltaY)), i - 1, nextJ, k, 2);
	maxwellEquationMatrix[i][j][k][1].push_back(element);

	//E i j+1 k+1
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i][nextJ][nextK].matrix[2][0] / (4 * deltaY * deltaZ)), i, nextJ, nextK, 0);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i][nextJ][nextK].matrix[2][1] / (4 * deltaY * deltaZ)), i, nextJ, nextK, 1);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i][nextJ][nextK].matrix[2][2] / (4 * deltaY * deltaZ)), i, nextJ, nextK, 2);
	maxwellEquationMatrix[i][j][k][1].push_back(element);

	//E i j-1 k-1
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i][prevJ][prevK].matrix[2][0] / (4 * deltaY * deltaZ)), i, prevJ, prevK, 0);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i][prevJ][prevK].matrix[2][1] / (4 * deltaY * deltaZ)), i, prevJ, prevK, 1);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i][prevJ][prevK].matrix[2][2] / (4 * deltaY * deltaZ)), i, prevJ, prevK, 2);
	maxwellEquationMatrix[i][j][k][1].push_back(element);

	//E i j+1 k-1
	element = MatrixElement(cthetadt2 * (dielectricTensor[i][nextJ][prevK].matrix[2][0] / (4 * deltaY * deltaZ)), i, nextJ, prevK, 0);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i][nextJ][prevK].matrix[2][1] / (4 * deltaY * deltaZ)), i, nextJ, prevK, 1);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i][nextJ][prevK].matrix[2][2] / (4 * deltaY * deltaZ)), i, nextJ, prevK, 2);
	maxwellEquationMatrix[i][j][k][1].push_back(element);

	//E i j-1 k+1
	element = MatrixElement(cthetadt2 * (dielectricTensor[i][prevJ][nextK].matrix[2][0] / (4 * deltaY * deltaZ)), i, prevJ, nextK, 0);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i][prevJ][nextK].matrix[2][1] / (4 * deltaY * deltaZ)), i, prevJ, nextK, 1);
	maxwellEquationMatrix[i][j][k][1].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i][prevJ][nextK].matrix[2][2] / (4 * deltaY * deltaZ)), i, prevJ, nextK, 2);
	maxwellEquationMatrix[i][j][k][1].push_back(element);

	double value = 0;
	for(int m = 0; m < maxwellEquationMatrix[i][j][k][1].size(); ++m) {
		value += maxwellEquationMatrix[i][j][k][1][m].value*Efield[maxwellEquationMatrix[i][j][k][1][m].i][maxwellEquationMatrix[i][j][k][1][m].j][maxwellEquationMatrix[i][j][k][1][m].k][maxwellEquationMatrix[i][j][k][1][m].l];
	}
	value = value - rightPart.y;
}

void Simulation::createInternalEquationZ(int i, int j, int k, Vector3d& rightPart) {
	int prevJ = j - 1;
	int prevK = k - 1;
	int nextJ = j + 1;
	int nextK = k + 1;
	if (prevJ < 0) {
		prevJ = ynumber - 1;
	}
	if (prevK < 0) {
		prevK = znumber - 1;
	}
	if (nextJ >= ynumber) {
		nextJ = 0;
	}
	if (nextK >= znumber) {
		nextK = 0;
	}
	double cthetadt2 = sqr(speed_of_light_normalized * theta * deltaT);

	MatrixElement element = MatrixElement(1 + dielectricTensor[i][j][k].matrix[2][2] + cthetadt2 * ((2 / (deltaX * deltaX)) + (2 / (deltaY * deltaY)) + (2 * (1 + 0.25 * dielectricTensor[i][j][k].matrix[2][2]) / (deltaZ * deltaZ))), i, j, k, 2);

	//E i j k
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(dielectricTensor[i][j][k].matrix[2][0] + cthetadt2 * ((0.5 * dielectricTensor[i][j][k].matrix[2][0] / (deltaZ * deltaZ))), i, j, k, 0);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(dielectricTensor[i][j][k].matrix[2][1] + cthetadt2 * ((0.5 * dielectricTensor[i][j][k].matrix[2][1] / (deltaZ * deltaZ))), i, j, k, 1);
	maxwellEquationMatrix[i][j][k][2].push_back(element);

	//E i j k+1
	element = MatrixElement(cthetadt2 * (-(0.5 * dielectricTensor[i][j][nextK].matrix[2][0]) / (deltaZ * deltaZ)), i, j, nextK, 0);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(cthetadt2 * (-(0.5 * dielectricTensor[i][j][nextK].matrix[2][1]) / (deltaZ * deltaZ)), i, j, nextK, 1);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(cthetadt2 * (-(1.0 + 0.5 * dielectricTensor[i][j][nextK].matrix[2][2]) / (deltaZ * deltaZ)), i, j, nextK, 2);
	maxwellEquationMatrix[i][j][k][2].push_back(element);

	//E i j k-1
	element = MatrixElement(cthetadt2 * (-(0.5 * dielectricTensor[i][j][prevK].matrix[2][0]) / (deltaZ * deltaZ)), i, j, prevK, 0);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(cthetadt2 * (-(0.5 * dielectricTensor[i][j][prevK].matrix[2][1]) / (deltaZ * deltaZ)), i, j, prevK, 1);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(cthetadt2 * (-(1.0 + 0.5 * dielectricTensor[i][j][prevK].matrix[2][2]) / (deltaZ * deltaZ)), i, j, prevK, 2);
	maxwellEquationMatrix[i][j][k][2].push_back(element);

	//Ez i+1 j k
	if (i < xnumber - 1) {
		element = MatrixElement(-cthetadt2 / (deltaX * deltaX), i + 1, j, k, 2);
		maxwellEquationMatrix[i][j][k][2].push_back(element);
	} else {
		rightPart.z += (cthetadt2 / (deltaX * deltaX))*E0.z;
	}

	//Ez i-1 j k
	element = MatrixElement(-cthetadt2 / (deltaX * deltaX), i - 1, j, k, 2);
	maxwellEquationMatrix[i][j][k][2].push_back(element);

	//Ez i j+1 k
	element = MatrixElement(-cthetadt2 / (deltaY * deltaY), i, nextJ, k, 2);
	maxwellEquationMatrix[i][j][k][2].push_back(element);

	//Ez i j+1 k
	element = MatrixElement(-cthetadt2 / (deltaY * deltaY), i, prevJ, k, 2);
	maxwellEquationMatrix[i][j][k][2].push_back(element);

	//E i+1 j k+1
	if (i < xnumber - 1) {
		element = MatrixElement(-cthetadt2 * (dielectricTensor[i + 1][j][nextK].matrix[0][0] / (4 * deltaX * deltaZ)), i + 1, j, nextK, 0);
		maxwellEquationMatrix[i][j][k][2].push_back(element);
		element = MatrixElement(-cthetadt2 * (dielectricTensor[i + 1][j][nextK].matrix[0][1] / (4 * deltaX * deltaZ)), i + 1, j, nextK, 1);
		maxwellEquationMatrix[i][j][k][2].push_back(element);
		element = MatrixElement(-cthetadt2 * (dielectricTensor[i + 1][j][nextK].matrix[0][2] / (4 * deltaX * deltaZ)), i + 1, j, nextK, 2);
		maxwellEquationMatrix[i][j][k][2].push_back(element);
	} else {
		rightPart.z += cthetadt2 * (dielectricTensor[i + 1][j][nextK].matrix[0][0] / (4 * deltaX * deltaZ))*E0.x;
		rightPart.z += cthetadt2 * (dielectricTensor[i + 1][j][nextK].matrix[0][1] / (4 * deltaX * deltaZ))*E0.y;
		rightPart.z += cthetadt2 * (dielectricTensor[i + 1][j][nextK].matrix[0][2] / (4 * deltaX * deltaZ))*E0.z;
	}

	//E i-1 j k-1
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i - 1][j][prevK].matrix[0][0] / (4 * deltaX * deltaZ)), i - 1, j, prevK, 0);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i - 1][j][prevK].matrix[0][1] / (4 * deltaX * deltaZ)), i - 1, j, prevK, 1);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i - 1][j][prevK].matrix[0][2] / (4 * deltaX * deltaZ)), i - 1, j, prevK, 2);
	maxwellEquationMatrix[i][j][k][2].push_back(element);

	//E i+1 j k-1
	if (i < xnumber - 1) {
		element = MatrixElement(cthetadt2 * (dielectricTensor[i + 1][j][prevK].matrix[0][0] / (4 * deltaX * deltaZ)), i + 1, j, prevK, 0);
		maxwellEquationMatrix[i][j][k][2].push_back(element);
		element = MatrixElement(cthetadt2 * (dielectricTensor[i + 1][j][prevK].matrix[0][1] / (4 * deltaX * deltaZ)), i + 1, j, prevK, 1);
		maxwellEquationMatrix[i][j][k][2].push_back(element);
		element = MatrixElement(cthetadt2 * (dielectricTensor[i + 1][j][prevK].matrix[0][2] / (4 * deltaX * deltaZ)), i + 1, j, prevK, 2);
		maxwellEquationMatrix[i][j][k][2].push_back(element);
	} else {
		rightPart.z -= cthetadt2 * (dielectricTensor[i + 1][j][prevK].matrix[0][0] / (4 * deltaX * deltaZ))*E0.x;
		rightPart.z -= cthetadt2 * (dielectricTensor[i + 1][j][prevK].matrix[0][1] / (4 * deltaX * deltaZ))*E0.y;
		rightPart.z -= cthetadt2 * (dielectricTensor[i + 1][j][prevK].matrix[0][2] / (4 * deltaX * deltaZ))*E0.z;
	}

	//E i-1 j k+1
	element = MatrixElement(cthetadt2 * (dielectricTensor[i - 1][j][nextK].matrix[0][0] / (4 * deltaX * deltaZ)), i - 1, j, nextK, 0);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i - 1][j][nextK].matrix[0][1] / (4 * deltaX * deltaZ)), i - 1, j, nextK, 1);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i - 1][j][nextK].matrix[0][2] / (4 * deltaX * deltaZ)), i - 1, j, nextK, 2);
	maxwellEquationMatrix[i][j][k][2].push_back(element);

	//E i j+1 k+1
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i][nextJ][nextK].matrix[1][0] / (4 * deltaY * deltaZ)), i, nextJ, nextK, 0);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i][nextJ][nextK].matrix[1][1] / (4 * deltaY * deltaZ)), i, nextJ, nextK, 1);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i][nextJ][nextK].matrix[1][2] / (4 * deltaY * deltaZ)), i, nextJ, nextK, 2);
	maxwellEquationMatrix[i][j][k][2].push_back(element);

	//E i j-1 k-1
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i][prevJ][prevK].matrix[1][0] / (4 * deltaY * deltaZ)), i, prevJ, prevK, 0);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i][prevJ][prevK].matrix[1][1] / (4 * deltaY * deltaZ)), i, prevJ, prevK, 1);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(-cthetadt2 * (dielectricTensor[i][prevJ][prevK].matrix[1][2] / (4 * deltaY * deltaZ)), i, prevJ, prevK, 2);
	maxwellEquationMatrix[i][j][k][2].push_back(element);

	//E i j+1 k-1
	element = MatrixElement(cthetadt2 * (dielectricTensor[i][nextJ][prevK].matrix[1][0] / (4 * deltaY * deltaZ)), i, nextJ, prevK, 0);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i][nextJ][prevK].matrix[1][1] / (4 * deltaY * deltaZ)), i, nextJ, prevK, 1);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i][nextJ][prevK].matrix[1][2] / (4 * deltaY * deltaZ)), i, nextJ, prevK, 2);
	maxwellEquationMatrix[i][j][k][2].push_back(element);

	//E i j-1 k+1
	element = MatrixElement(cthetadt2 * (dielectricTensor[i][prevJ][nextK].matrix[1][0] / (4 * deltaY * deltaZ)), i, prevJ, nextK, 0);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i][prevJ][nextK].matrix[1][1] / (4 * deltaY * deltaZ)), i, prevJ, nextK, 1);
	maxwellEquationMatrix[i][j][k][2].push_back(element);
	element = MatrixElement(cthetadt2 * (dielectricTensor[i][prevJ][nextK].matrix[1][2] / (4 * deltaY * deltaZ)), i, prevJ, nextK, 2);
	maxwellEquationMatrix[i][j][k][2].push_back(element);

	double value = 0;
	for(int m = 0; m < maxwellEquationMatrix[i][j][k][2].size(); ++m) {
		value += maxwellEquationMatrix[i][j][k][2][m].value*Efield[maxwellEquationMatrix[i][j][k][2][m].i][maxwellEquationMatrix[i][j][k][2][m].j][maxwellEquationMatrix[i][j][k][2][m].k][maxwellEquationMatrix[i][j][k][2][m].l];
	}
	value = value - rightPart.z;
}

void Simulation::evaluateMagneticField() {
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				Vector3d rotE = evaluateRotE(i, j, k);
				newBfield[i][j][k] = Bfield[i][j][k] + rotE * speed_of_light_normalized * deltaT;
			}
		}
	}
}

void Simulation::updateBoundaries() {
	for (int i = 0; i < xnumber; ++i) {
		//periodic by Y
		for (int k = 0; k < znumber; ++k) {
			tempEfield[i][ynumber][k] = tempEfield[i][0][k];
		}
		//periodic by Z
		//note j <= number because corner point

		for (int j = 0; j <= ynumber; ++j) {
			tempEfield[i][j][znumber] = tempEfield[i][j][0];
		}
	}

	for (int j = 0; j <= ynumber; ++j) {
		for (int k = 0; k <= znumber; ++k) {
			tempEfield[xnumber][j][k] = E0;
		}
	}
}

void Simulation::updateBoundariesOldField() {
	for (int i = 0; i < xnumber; ++i) {
		//periodic by Y
		for (int k = 0; k < znumber; ++k) {
			Efield[i][ynumber][k] = Efield[i][0][k];
		}
		//periodic by Z
		//note j <= number because corner point

		for (int j = 0; j <= ynumber; ++j) {
			Efield[i][j][znumber] = Efield[i][j][0];
		}
	}

	for (int j = 0; j <= ynumber; ++j) {
		for (int k = 0; k <= znumber; ++k) {
			Efield[xnumber][j][k] = E0;
		}
	}
}

double Simulation::evaluateDivFlux(int i, int j, int k) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			exit(0);
		}

		if (i >= xnumber) {
			printf("x >= xnumber\n");
			exit(0);
		}

		if (j < 0) {
			printf("j < 0\n");
			exit(0);
		}

		if (j >= ynumber) {
			printf("j >= ynumber\n");
			exit(0);
		}

		if (k < 0) {
			printf("k < 0\n");
			exit(0);
		}

		if (k >= znumber) {
			printf("k >= znumber\n");
			exit(0);
		}
	}


	double rfluxX = (electricFlux[i + 1][j][k].x + electricFlux[i + 1][j + 1][k].x + electricFlux[i + 1][j][k + 1].x + electricFlux[i + 1][j + 1][k + 1].x) / 4;
	double lfluxX = (electricFlux[i][j][k].x + electricFlux[i][j + 1][k].x + electricFlux[i][j][k + 1].x + electricFlux[i][j + 1][k + 1].x) / 4;

	double rfluxY = (electricFlux[i][j + 1][k].y + electricFlux[i + 1][j + 1][k].y + electricFlux[i][j + 1][k + 1].y + electricFlux[i + 1][j + 1][k + 1].y) / 4;
	double lfluxY = (electricFlux[i][j][k].y + electricFlux[i + 1][j][k].y + electricFlux[i][j][k + 1].y + electricFlux[i + 1][j][k + 1].y) / 4;

	double rfluxZ = (electricFlux[i][j][k + 1].z + electricFlux[i + 1][j][k + 1].z + electricFlux[i][j + 1][k + 1].z + electricFlux[i + 1][j + 1][k + 1].z) / 4;
	double lfluxZ = (electricFlux[i][j][k].z + electricFlux[i + 1][j][k].z + electricFlux[i][j + 1][k].z + electricFlux[i + 1][j + 1][k].z) / 4;

	return ((rfluxX - lfluxX) / deltaX) + ((rfluxY - lfluxY) / deltaY) + ((rfluxZ - lfluxZ) / deltaZ);
}

Vector3d Simulation::evaluateRotB(int i, int j, int k) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			exit(0);
		}

		if (i > xnumber) {
			printf("x > xnumber\n");
			exit(0);
		}

		if (j < 0) {
			printf("j < 0\n");
			exit(0);
		}

		if (j > ynumber) {
			printf("j > ynumber\n");
			exit(0);
		}

		if (k < 0) {
			printf("k < 0\n");
			exit(0);
		}

		if (k > znumber) {
			printf("k > znumber\n");
			exit(0);
		}
	}

	Vector3d BrightX;
	Vector3d BleftX;
	Vector3d BrightY;
	Vector3d BleftY;
	Vector3d BrightZ;
	Vector3d BleftZ;

	if (i == 0) {
	} else {
		BrightX = (getBfield(i, j, k) + getBfield(i, j - 1, k) + getBfield(i, j - 1, k - 1) + getBfield(i, j - 1, k - 1)) / 4;
		BleftX = (getBfield(i - 1, j, k) + getBfield(i - 1, j - 1, k) + getBfield(i - 1, j - 1, k - 1) + getBfield(i - 1, j - 1, k - 1)) / 4;

		BrightY = (getBfield(i, j, k) + getBfield(i - 1, j, k) + getBfield(i, j, k - 1) + getBfield(i - 1, j, k - 1)) / 4;
		BleftY = (getBfield(i, j - 1, k) + getBfield(i - 1, j - 1, k) + getBfield(i, j - 1, k - 1) + getBfield(i - 1, j - 1, k - 1)) / 4;

		BrightZ = (getBfield(i, j, k) + getBfield(i - 1, j, k) + getBfield(i, j - 1, k) + getBfield(i - 1, j - 1, k)) / 4;
		BleftZ = (getBfield(i, j, k - 1) + getBfield(i - 1, j, k - 1) + getBfield(i, j - 1, k - 1) + getBfield(i - 1, j - 1, k - 1)) / 4;
	}


	double x = 0;
	double y = 0;
	double z = 0;

	x = ((BrightY.z - BleftY.z) / deltaY) - ((BrightZ.y - BleftZ.y) / deltaZ);
	y = ((BrightZ.x - BleftZ.x) / deltaZ) - ((BrightX.z - BleftX.z) / deltaX);
	z = ((BrightX.y - BleftX.y) / deltaX) - ((BrightY.x - BleftY.x) / deltaY);

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateRotE(int i, int j, int k) {
	if (debugMode) {
		if (i < 0) {
			printf("i < 0\n");
			exit(0);
		}

		if (i >= xnumber) {
			printf("x >= xnumber\n");
			exit(0);
		}

		if (j < 0) {
			printf("j < 0\n");
			exit(0);
		}

		if (j >= ynumber) {
			printf("j >= ynumber\n");
			exit(0);
		}

		if (k < 0) {
			printf("k < 0\n");
			exit(0);
		}

		if (k >= znumber) {
			printf("k >= znumber\n");
			exit(0);
		}
	}

	double x = 0;
	double y = 0;
	double z = 0;

	Vector3d ErightX = (tempEfield[i + 1][j][k] + tempEfield[i + 1][j + 1][k] + tempEfield[i + 1][j][k + 1] + tempEfield[i + 1][j + 1][k + 1]) / 4;
	Vector3d EleftX = (tempEfield[i][j][k] + tempEfield[i][j + 1][k] + tempEfield[i][j][k + 1] + tempEfield[i][j + 1][k + 1]) / 4;

	Vector3d ErightY = (tempEfield[i][j + 1][k] + tempEfield[i + 1][j + 1][k] + tempEfield[i][j + 1][k + 1] + tempEfield[i + 1][j + 1][k + 1]) / 4;
	Vector3d EleftY = (tempEfield[i][j][k] + tempEfield[i + 1][j][k] + tempEfield[i][j][k + 1] + tempEfield[i + 1][j][k + 1]) / 4;

	Vector3d ErightZ = (tempEfield[i][j][k + 1] + tempEfield[i + 1][j][k + 1] + tempEfield[i][j + 1][k + 1] + tempEfield[i + 1][j + 1][k + 1]) / 4;
	Vector3d EleftZ = (tempEfield[i][j][k] + tempEfield[i + 1][j][k] + tempEfield[i][j + 1][k] + tempEfield[i + 1][j + 1][k]) / 4;

	x = ((ErightY.z - EleftY.z) / deltaY) - ((ErightZ.y - EleftZ.y) / deltaZ);
	y = ((ErightZ.x - EleftZ.x) / deltaZ) - ((ErightX.z - EleftX.z) / deltaX);
	z = ((ErightX.y - EleftX.y) / deltaX) - ((ErightY.x - EleftY.x) / deltaY);

	return Vector3d(x, y, z);
}

Vector3d Simulation::evaluateDivPressureTensor(int i, int j, int k) {
	Vector3d result = Vector3d(0, 0, 0);

	Matrix3d tensorDerX = (getPressureTensor(i, j, k) + getPressureTensor(i, j - 1, k) + getPressureTensor(i, j, k - 1) + getPressureTensor(i, j - 1, k - 1) - getPressureTensor(i - 1, j, k) - getPressureTensor(i - 1, j - 1, k) - getPressureTensor(i - 1, j, k - 1) - getPressureTensor(i - 1, j - 1, k - 1)) / (4 * deltaX);
	Matrix3d tensorDerY = (getPressureTensor(i, j, k) + getPressureTensor(i - 1, j, k) + getPressureTensor(i, j, k - 1) + getPressureTensor(i - 1, j, k - 1) - getPressureTensor(i, j - 1, k) - getPressureTensor(i - 1, j - 1, k) - getPressureTensor(i, j - 1, k - 1) - getPressureTensor(i - 1, j - 1, k - 1)) / (4 * deltaY);
	Matrix3d tensorDerZ = (getPressureTensor(i, j, k) + getPressureTensor(i - 1, j, k) + getPressureTensor(i, j - 1, k) + getPressureTensor(i - 1, j - 1, k) - getPressureTensor(i, j, k - 1) - getPressureTensor(i - 1, j, k - 1) - getPressureTensor(i, j - 1, k - 1) - getPressureTensor(i - 1, j - 1, k - 1)) / (4 * deltaZ);

	result.x = tensorDerX.matrix[0][0] + tensorDerY.matrix[1][0] + tensorDerZ.matrix[2][0];
	result.y = tensorDerX.matrix[0][1] + tensorDerY.matrix[1][1] + tensorDerZ.matrix[2][1];
	result.z = tensorDerX.matrix[0][2] + tensorDerY.matrix[1][2] + tensorDerZ.matrix[2][2];

	return result;
}

Vector3d Simulation::evaluateGradDensity(int i, int j, int k) {

	double densityRightX = (getDensity(i + 1, j, k) + getDensity(i + 1, j + 1, k) + getDensity(i + 1, j, k + 1) + getDensity(i + 1, j + 1, k + 1)) / 4;
	double densityLeftX = (getDensity(i, j, k) + getDensity(i, j + 1, k) + getDensity(i, j, k + 1) + getDensity(i, j + 1, k + 1)) / 4;

	double densityRightY = (getDensity(i, j + 1, k) + getDensity(i + 1, j + 1, k) + getDensity(1, j + 1, k + 1) + getDensity(i + 1, j + 1, k + 1)) / 4;
	double densityLeftY = (getDensity(i, j, k) + getDensity(i + 1, j, k) + getDensity(1, j, k + 1) + getDensity(i + 1, j, k + 1)) / 4;

	double densityRightZ = (getDensity(i, j, k + 1) + getDensity(i + 1, j, k + 1) + getDensity(i, j + 1, k + 1) + getDensity(i + 1, j + 1, k + 1)) / 4;
	double densityLeftZ = (getDensity(i, j, k) + getDensity(i + 1, j, k) + getDensity(i, j + 1, k) + getDensity(i + 1, j + 1, k)) / 4;

	double x = (densityRightX - densityLeftX) / deltaX;
	double y = (densityRightY - densityLeftY) / deltaY;
	double z = (densityRightZ - densityLeftZ) / deltaZ;

	return Vector3d(x, y, z);
}