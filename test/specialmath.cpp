#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "specialmath.h"


double** arnoldiIterations(double** matrix, double** outHessenbergMatrix, int n, double** prevBasis, double** prevHessenbergMatrix){
	double** resultBasis = new double*[n];
	for(int m = 0; m < n-1; ++m){
		resultBasis[m] = prevBasis[m];
	}
	delete[] prevBasis;

	for(int i = 0; i < n; ++i){
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

	double* tempVector = multiplyMatrixVector(matrix, resultBasis[n-2]);

	for(int m = 0; m < n-1; ++m){
		outHessenbergMatrix[m][n-2] = scalarMultiplyLargeVectors(resultBasis[m], tempVector);
		for(int i = 0; i < number ; ++i){
			tempVector[i] -= outHessenbergMatrix[m][n-1]*resultBasis[m][i];
		}
	}
	outHessenbergMatrix[n-1][n-2] = sqrt(scalarMultiplyLargeVectors(tempVector, tempVector));
	for(int i = 0; i < number; ++i){
		tempVector[i] /= outHessenbergMatrix[n-1][n-2];
	}

	resultBasis[n-1] = tempVector;

	return resultBasis;
}

void generalizedMinimalResidualMethod(double** matrix, double* rightPart){
	double norm = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart));

	for(int i = 0; i < number; ++i){
		rightPart[i] /= norm;
		for(int m = 0; m < number; ++m){
			matrix[i][m] /= norm;
		}
	}

	double maxError = 0.0001;

	double** hessenbergMatrix;
	double** newHessenbergMatrix = NULL;
	hessenbergMatrix = new double*[1];
	hessenbergMatrix[0] = new double[1];

	double** Qmatrix = new double*[2];
	double** Rmatrix = new double*[2];
	double** oldQmatrix = new double*[2];
	double** oldRmatrix = new double*[2];

	for(int i = 0; i < 2; ++i){
		Qmatrix[i] = new double[2];
		oldQmatrix[i] = new double[2];
	}

	Rmatrix[0] = new double[1];
	Rmatrix[1] = new double[1];
	oldRmatrix[0] = new double[1];
	oldRmatrix[1] = new double[1];

	double** basis = new double*[1];
	basis[0] = new double[number];
	for(int i = 0; i < number; ++i){
		basis[0][i] = rightPart[i];
	}
	double** newBasis;

	int n = 2;
	double beta = 1.0;
	double error = beta;
	double* y = new double[1];

	double rho;
	double sigma;
	double cosn;
	double sinn;
	double module;

	while(error > maxError && n < number){
		printf("iteration %d\n", n);
		newHessenbergMatrix = new double*[n];
		for(int i = 0; i < n; ++i){
			newHessenbergMatrix[i] = new double[n-1];
		}
		newBasis = arnoldiIterations(matrix, newHessenbergMatrix, n, basis, hessenbergMatrix);

		hessenbergMatrix = newHessenbergMatrix;
		basis = newBasis;

		if(n == 2){
			rho = hessenbergMatrix[0][0];
			sigma = hessenbergMatrix[1][0];

			module = sqrt(rho*rho + sigma*sigma);
			
			cosn = rho/module;
			sinn = sigma/module;

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
			for(int i = 0; i < n; ++i){
				Rmatrix[i] = new double[n-1];
				if(i < n-1){
					for(int j = 0; j < n-2; ++j){
						Rmatrix[i][j] = oldRmatrix[i][j];
					}
				} else {
					for(int j = 0; j < n-2; ++j){
						Rmatrix[i][j] = 0;
					}
				}
			}

			Qmatrix = new double*[n];
			for(int i = 0; i < n;++i){
				Qmatrix[i] = new double[n];
				if(i < n-1){
					for(int j = 0; j < n-1; ++j){
						Qmatrix[i][j] = oldQmatrix[i][j];
					}
					Qmatrix[i][n-1] = 0;
				} else {
					for(int j = 0; j < n-1; ++j){
						Qmatrix[i][j] = 0;
					}
					Qmatrix[n-1][n-1] = 1;
				}
			}

			for(int i = 0; i < n; ++i){
				Rmatrix[i][n-2] = 0;
				for(int j = 0; j < n; ++j){
					Rmatrix[i][n-2] += Qmatrix[i][j]*hessenbergMatrix[j][n-2];
				}
			}
			rho = Rmatrix[n-2][n-2];
			sigma = Rmatrix[n-1][n-2];

			module = sqrt(rho*rho + sigma*sigma);

			cosn = rho/module;
			sinn = sigma/module;

			Rmatrix[n-2][n-2] = module;
			Rmatrix[n-1][n-2] = 0;

			for(int j = 0; j < n-1; ++j){
				Qmatrix[n-2][j] = cosn*oldQmatrix[n-2][j];
				Qmatrix[n-1][j] = - -sinn*oldQmatrix[n-2][j];
			}
			Qmatrix[n-2][n-1] = sinn;
			Qmatrix[n-1][n-1] = cosn;
		}

		delete[] y;
		y = new double[n-1];

		for(int i = n-2; i >=0; --i){
			y[i] = beta*Qmatrix[i][0];
			for(int j = n-2; j > i; --j){
				y[i] -= Rmatrix[i][j]*y[j];
			}
			y[i] /= Rmatrix[i][i];
		}

		error = abs(beta*Qmatrix[n-1][0]);

		for(int i = 0; i < n-1; ++i){
			delete[] oldQmatrix[i];
			delete[] oldRmatrix[i];
		}
		delete[] oldQmatrix;
		delete[] oldRmatrix;

		oldQmatrix = Qmatrix;
		oldRmatrix = Rmatrix;

		n++;
	}

	n = n-1;

	//out result

	for(int i = 0; i < number; ++i){
		double a = 0;
		for(int m = 0; m < n; ++m){
			a += basis[m][i]*y[m];
		}
		printf("%lf\n", a);
	}

	for(int i = 0; i < n; ++i){
		delete[] Qmatrix[i];
		delete[] Rmatrix[i];
		delete[] hessenbergMatrix[i];
	}
	delete[] Qmatrix;
	delete[] Rmatrix;
	delete[] hessenbergMatrix;

	for(int m = 0; m < n; ++m){
		delete[] basis[m];
	}
	delete[] basis;

	delete[] y;
}

double scalarMultiplyLargeVectors(double* a, double* b){
	double result = 0;
	for(int i = 0; i < number; ++i){
		result += a[i]*b[i];
	}
	return result;
}

double* multiplyMatrixVector(double** matrix, double* vector){
	double* result = new double[number];
	for(int i = 0; i < number; ++i){
		result[i] = 0;
		for(int j = 0; j < number; ++j){
			result[i] += matrix[i][j]*vector[j];
		}
	}

	return result;
}