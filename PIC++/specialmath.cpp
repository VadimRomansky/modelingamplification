#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "simulation.h"
#include "util.h"
#include "constants.h"

double Simulation::evaluateError(double** hessenbergMatrix, double* vector, double beta, int n) {
	double* resVector = new double[n + 1];

	for (int i = 0; i < n + 1; ++i) {
		resVector[i] = 0;
		for (int j = 0; j < n; ++j) {
			resVector[i] += hessenbergMatrix[i][j] * vector[j];
		}
		if (i == 0) {
			resVector[i] -= beta;
		}
	}

	double norm = 0;
	for (int i = 0; i < n + 1; ++i) {
		norm += resVector[i] * resVector[i];
	}

	delete[] resVector;

	return sqrt(norm);
}

double**** Simulation::multiplySpecialMatrixVector(double**** vector) {
	double**** result = new double***[xnumber];
	for (int i = 0; i < xnumber; ++i) {
		result[i] = new double**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			result[i][j] = new double*[znumber];
			for (int k = 0; k < znumber; ++k) {
				result[i][j][k] = new double[3];
				for (int l = 0; l < 3; ++l) {
					result[i][j][k][l] = 0;
					for (int m = 0; m < maxwellEquationMatrix[i][j][k][l].size(); ++i) {
						MatrixElement element = maxwellEquationMatrix[i][j][k][l][m];

						result[i][j][k][l] += element.value * vector[element.i][element.j][element.k][element.l];
					}
				}
			}
		}
	}

	return result;
}

double***** Simulation::arnoldiIterations(double** outHessenbergMatrix, int n, double***** prevBasis, double** prevHessenbergMatrix) {
	double***** resultBasis = new double****[n];
	for (int m = 0; m < n - 1; ++m) {
		resultBasis[m] = prevBasis[m];
	}
	delete[] prevBasis;

	for (int i = 0; i < n; ++i) {
		if (i < n - 1) {
			for (int j = 0; j < n - 2; ++j) {
				outHessenbergMatrix[i][j] = prevHessenbergMatrix[i][j];
			}
			delete[] prevHessenbergMatrix[i];
		} else {
			for (int j = 0; j < n - 2; ++j) {
				outHessenbergMatrix[i][j] = 0;
			}
		}
	}
	delete[] prevHessenbergMatrix;

	double**** tempVector = multiplySpecialMatrixVector(resultBasis[n - 2]);

	for (int m = 0; m < n - 1; ++m) {
		outHessenbergMatrix[m][n - 2] = scalarMultiplyLargeVectors(resultBasis[m], tempVector);
		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber + 1; ++j) {
				for (int k = 0; k < znumber + 1; ++k) {
					for (int l = 0; l < 3; ++l) {
						tempVector[i][j][k][l] -= outHessenbergMatrix[m][n - 2] * resultBasis[m][i][j][k][l];
					}
				}
			}
		}
	}
	outHessenbergMatrix[n - 1][n - 2] = sqrt(scalarMultiplyLargeVectors(tempVector, tempVector));
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					tempVector[i][j][k][l] /= outHessenbergMatrix[n - 1][n - 2];
				}
			}
		}
	}

	resultBasis[n - 1] = tempVector;

	return resultBasis;
}

void Simulation::generalizedMinimalResidualMethod() {
	double norm = sqrt(scalarMultiplyLargeVectors(maxwellEquationRightPart, maxwellEquationRightPart));

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				for (int l = 0; l < 3; ++l) {
					maxwellEquationRightPart[i][j][k][l] /= norm;
					for (int m = 0; m < maxwellEquationMatrix[i][j][k][l].size(); ++m) {
						maxwellEquationMatrix[i][j][k][l][m].value /= norm;
					}
				}
			}
		}
	}

	double maxError = E0.getNorm() / 100;

	double** hessenbergMatrix;
	double** newHessenbergMatrix;
	hessenbergMatrix = new double*[1];
	hessenbergMatrix[0] = new double[1];

	double** Qmatrix = new double*[2];
	double** Rmatrix = new double*[2];
	double** oldQmatrix = new double*[2];
	double** oldRmatrix = new double*[2];

	for (int i = 0; i < 2; ++i) {
		Qmatrix[i] = new double[2];
		oldQmatrix[i] = new double[2];
	}

	Rmatrix[0] = new double[1];
	oldRmatrix[0] = new double[1];

	double***** basis = new double****[1];
	basis[0] = new double***[xnumber + 1];
	for (int i = 0; i < xnumber + 1; ++i) {
		basis[0][i] = new double**[xnumber + 1];
		for (int j = 0; j < ynumber + 1; ++j) {
			basis[0][i][j] = new double*[ynumber + 1];
			for (int k = 0; k < znumber + 1; ++k) {
				basis[0][i][j][k] = new double[znumber + 1];
				for (int l = 0; l < 3; ++l) {
					basis[0][i][j][k][l] = maxwellEquationRightPart[i][j][k][l];
				}
			}
		}
	}
	double***** newBasis;

	int n = 2;
	double beta = 1.0;
	double error = beta;
	double* y = new double[1];

	double rho;
	double sigma;
	double cosn;
	double sinn;
	double module;

	while (error > maxError && n < maxGMRESIterations) {
		printf("GMRES iteration %d\n", n);
		newHessenbergMatrix = new double*[n];
		for (int i = 0; i < n; ++i) {
			newHessenbergMatrix[i] = new double[n - 1];
		}
		newBasis = arnoldiIterations(newHessenbergMatrix, n, basis, hessenbergMatrix);

		hessenbergMatrix = newHessenbergMatrix;
		basis = newBasis;

		if (n == 2) {
			rho = hessenbergMatrix[0][0];
			sigma = hessenbergMatrix[1][0];

			module = sqrt(rho * rho + sigma * sigma);

			cosn = rho / module;
			sinn = sigma / module;

			Qmatrix[0][0] = cosn;
			Qmatrix[0][1] = sinn;
			Qmatrix[1][0] = -sinn;
			Qmatrix[1][1] = cosn;

			oldQmatrix[0][0] = Qmatrix[0][0];
			oldQmatrix[0][1] = Qmatrix[0][1];
			oldQmatrix[1][0] = Qmatrix[1][0];
			oldQmatrix[1][1] = Qmatrix[1][1];

			Rmatrix[0][0] = module;
			Rmatrix[1][0] = 0;

			oldRmatrix[0][0] = Rmatrix[0][0];
			oldRmatrix[1][0] = Rmatrix[1][0];

		} else {
			Rmatrix = new double*[n];
			for (int i = 0; i < n; ++i) {
				Rmatrix[i] = new double[n - 1];
				if (i < n - 1) {
					for (int j = 0; j < n - 2; ++j) {
						Rmatrix[i][j] = oldRmatrix[i][j];
					}
				} else {
					for (int j = 0; j < n - 2; ++j) {
						Rmatrix[i][j] = 0;
					}
				}
			}

			Qmatrix = new double*[n];
			for (int i = 0; i < n;++i) {
				Qmatrix[i] = new double[n];
				if (i < n - 1) {
					for (int j = 0; j < n - 1; ++j) {
						Qmatrix[i][j] = oldQmatrix[i][j];
					}
					Qmatrix[i][n - 1] = 0;
				} else {
					for (int j = 0; j < n - 1; ++j) {
						Qmatrix[i][j] = 0;
					}
					Qmatrix[n - 1][n - 1] = 1;
				}
			}

			for (int i = 0; i < n; ++i) {
				Rmatrix[i][n - 2] = 0;
				for (int j = 0; j < n; ++j) {
					Rmatrix[i][n - 2] += Qmatrix[i][j] * hessenbergMatrix[j][n - 2];
				}
			}
			rho = Rmatrix[n - 2][n - 2];
			sigma = Rmatrix[n - 1][n - 2];

			module = sqrt(rho * rho + sigma * sigma);

			cosn = rho / module;
			sinn = sigma / module;

			Rmatrix[n - 2][n - 2] = module;
			Rmatrix[n - 1][n - 2] = 0;

			for (int j = 0; j < n - 1; ++j) {
				Qmatrix[n - 2][j] = cosn * oldQmatrix[n - 2][j];
				Qmatrix[n - 1][j] = -sinn * oldQmatrix[n - 2][j];
			}
			Qmatrix[n - 2][n - 1] = sinn;
			Qmatrix[n - 1][n - 1] = cosn;
		}

		delete[] y;
		y = new double[n - 1];

		for (int i = n - 2; i >= 0; --i) {
			y[i] = beta * Qmatrix[i][0];
			for (int j = n - 2; j > i; --j) {
				y[i] -= Rmatrix[i][j] * y[j];
			}
			y[i] /= Rmatrix[i][i];
		}

		error = fabs(beta * Qmatrix[n - 1][0]);

		for (int i = 0; i < n - 1; ++i) {
			delete[] oldQmatrix[i];
			delete[] oldRmatrix[i];
		}
		delete[] oldQmatrix;
		delete[] oldRmatrix;

		oldQmatrix = Qmatrix;
		oldRmatrix = Rmatrix;

		n++;
	}

	n = n - 1;

	//out result

	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				tempEfield[i][j][k].x = 0;
				tempEfield[i][j][k].y = 0;
				tempEfield[i][j][k].z = 0;
				for (int m = 0; m < n; ++m) {
					tempEfield[i][j][k].x += basis[m][i][j][k][0] * y[m];
					tempEfield[i][j][k].y += basis[m][i][j][k][1] * y[m];
					tempEfield[i][j][k].z += basis[m][i][j][k][2] * y[m];
				}
			}
		}
	}

	for (int i = 0; i < n; ++i) {
		delete[] Qmatrix[i];
		delete[] Rmatrix[i];
		delete[] hessenbergMatrix[i];
	}
	delete[] Qmatrix;
	delete[] Rmatrix;
	delete[] hessenbergMatrix;

	for (int m = 0; m < n; ++m) {
		for (int i = 0; i < xnumber + 1; ++i) {
			for (int j = 0; j < ynumber + 1; ++j) {
				for (int k = 0; k < znumber + 1; ++k) {
					delete[] basis[m][i][j][k];
				}
				delete[] basis[m][i][j];
			}
			delete[] basis[m][i];
		}
		delete[] basis[m];
	}
	delete[] basis;

	delete[] y;
}

double Simulation::scalarMultiplyLargeVectors(double**** a, double**** b) {
	double result = 0;
	for (int i = 0; i < xnumber + 1; ++i) {
		for (int j = 0; j < ynumber + 1; ++j) {
			for (int k = 0; k < znumber + 1; ++k) {
				for (int l = 0; l < 3; ++i) {
					result += a[i][j][k][l] * b[i][j][k][l];
				}
			}
		}
	}
	return result;
}