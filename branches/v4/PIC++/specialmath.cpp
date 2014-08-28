#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "simulation.h"
#include "util.h"
#include "constants.h"

double Simulation::evaluateError(double** hessenbergMatrix, double* vector, double beta, int n){
	double* resVector = new double[n+1];

	for(int i = 0; i < n+1; ++i){
		resVector[i] = 0;
		for(int j = 0; j < n; ++j){
			resVector[i] += hessenbergMatrix[i][j]*vector[j];
		}
		if(i == 0) {
			resVector[i] -= beta;
		}
	}

	double norm = 0;
	for(int i = 0; i < n+1; ++i){
		norm += resVector[i]*resVector[i];
	}

	delete[] resVector;

	return sqrt(norm);
}

double**** Simulation::multiplySpecialMatrixVector(double**** vector){
}

double***** Simulation::arnoldiIterations(double** outHessenbergMatrix, int n, double***** prevBasis, double** prevHessenbergMatrix){
	double***** resultBasis = new double****[n];
	for(int m = 0; m < n-1; ++m){
		resultBasis[m] = prevBasis[m];
	}
	delete[] prevBasis;

	outHessenbergMatrix = new double*[n];
	for(int i = 0; i < n; ++i){
		outHessenbergMatrix[i] = new double[n-1];
		if(i < n - 1){
			for(int j = 0; j < n-2; ++j){
				outHessenbergMatrix[i][j] = prevHessenbergMatrix[i][j];
			}
			delete[] prevHessenbergMatrix[i];
		} else {
			for(int j= 0; j < n-2; ++j){
				outHessenbergMatrix[i][j] = 0;
			}
		}
	}
	delete[] prevHessenbergMatrix;

	double**** tempVector = multiplySpecialMatrixVector(resultBasis[n-2]);

	for(int m = 0; m < n; ++m){
		outHessenbergMatrix[m][n-2] = scalarMultiplyLargeVectors(resultBasis[m], tempVector);
		for(int i = 0; i < xnumber + 1; ++i){
			for(int j = 0; j < ynumber + 1; ++j){
				for(int k = 0; k < znumber + 1; ++k){
					for(int l = 0; l < 3; ++l){
						tempVector[i][j][k][l] -= outHessenbergMatrix[m][n-1]*resultBasis[m][i][j][k][l];
					}
				}
			}
		}
	}
	outHessenbergMatrix[n-1][n-2] = sqrt(scalarMultiplyLargeVectors(tempVector, tempVector));
	for(int i = 0; i < xnumber + 1; ++i){
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				for(int l = 0; l < 3; ++l){
					tempVector[i][j][k][l] /= outHessenbergMatrix[n-1][n-2];
				}
			}
		}
	}

	resultBasis[n-1] = tempVector;

	return resultBasis;
}

void Simulation::generalizedMinimalResidualMethod(){
	double norm = sqrt(scalarMultiplyLargeVectors(maxwellEquationRightPart, maxwellEquationRightPart));

	for(int i = 0; i < xnumber + 1; ++i){
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				for(int l = 0; l < 3; ++l){
					maxwellEquationRightPart[i][j][k][l] /= norm;
					for(int m = 0; m < 29; ++m){
						maxwellEquationMatrix[i][j][k][l][m] /= norm;
					}
				}
			}
		}
	}

	double maxError = E0.getNorm()/100;

	double** hessenbergMatrix;
	double** newHessenbergMatrix;
	hessenbergMatrix = new double*[1];
	hessenbergMatrix[0] = new double[1];

	double***** basis = new double****[1];
	basis[0] = new double***[xnumber + 1];
	for(int i = 0; i < xnumber + 1; ++i){
		basis[0][i] = new double**[xnumber + 1];
		for(int j = 0; j < ynumber + 1; ++j){
			basis[0][i][j] = new double*[ynumber + 1];
			for(int k = 0; k < znumber + 1; ++k){
				basis[0][i][j][k] = new double[znumber + 1];
				for(int l = 0; l < 3; ++l){
					basis[0][i][j][k][l] = maxwellEquationRightPart[i][j][k][l];
				}
			}
		}
	}
	double***** newBasis;

	int n = 2;
	double beta = 1.0;
	double error = beta;
	double* y;

	while(error > maxError && n < maxGMRESIterations){
		newBasis = arnoldiIterations(newHessenbergMatrix, n, basis, hessenbergMatrix);
		y = linearLeastSquares(newHessenbergMatrix, n);
		error = evaluateError(newHessenbergMatrix, y, beta, n);
		hessenbergMatrix = newHessenbergMatrix;
		basis = newBasis;
		n++;
	}

	n = n-1;
}

double Simulation::scalarMultiplyLargeVectors(double**** a, double**** b){
	double result = 0;
	for(int i = 0; i < xnumber + 1; ++i){
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				for(int l = 0; l < 3; ++i){
					result += a[i][j][k][l]*b[i][j][k][l];
				}
			}
		}
	}
	return result;
}