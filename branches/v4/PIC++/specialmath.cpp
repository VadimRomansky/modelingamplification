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

double**** Simulation::multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, Vector3d*** vector) {
	double**** result = new double***[xnumber];
	for (int i = 0; i < xnumber; ++i) {
		result[i] = new double**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			result[i][j] = new double*[znumber];
			for (int k = 0; k < znumber; ++k) {
				result[i][j][k] = new double[3];
				for (int l = 0; l < 3; ++l) {
					result[i][j][k][l] = 0;
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						MatrixElement element = matrix[i][j][k][l][m];

						result[i][j][k][l] += element.value * vector[element.i][element.j][element.k][element.l];
					}
				}
			}
		}
	}

	return result;
}

double**** Simulation::multiplySpecialMatrixVector(std::vector<MatrixElement>**** matrix, double**** vector) {
	double**** result = new double***[xnumber];
	for (int i = 0; i < xnumber; ++i) {
		result[i] = new double**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			result[i][j] = new double*[znumber];
			for (int k = 0; k < znumber; ++k) {
				result[i][j][k] = new double[3];
				for (int l = 0; l < 3; ++l) {
					result[i][j][k][l] = 0;
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						MatrixElement element = matrix[i][j][k][l][m];

						result[i][j][k][l] += element.value * vector[element.i][element.j][element.k][element.l];
					}
				}
			}
		}
	}

	return result;
}

double***** Simulation::arnoldiIterations(std::vector<MatrixElement>**** matrix,double** outHessenbergMatrix, int n, double***** prevBasis, double** prevHessenbergMatrix) {
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

	double**** tempVector = multiplySpecialMatrixVector(matrix, resultBasis[n - 2]);

	for (int m = 0; m < n - 1; ++m) {
		outHessenbergMatrix[m][n - 2] = scalarMultiplyLargeVectors(resultBasis[m], tempVector);
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < 3; ++l) {
						tempVector[i][j][k][l] -= outHessenbergMatrix[m][n - 2] * resultBasis[m][i][j][k][l];
						alertNaNOrInfinity(tempVector[i][j][k][l], "tempVector = NaN\n");
					}
				}
			}
		}
	}
	double a = scalarMultiplyLargeVectors(resultBasis[0], tempVector);
	outHessenbergMatrix[n - 1][n - 2] = sqrt(scalarMultiplyLargeVectors(tempVector, tempVector));
	if (outHessenbergMatrix[n - 1][n - 2] > 0) {
		for (int i = 0; i < xnumber ; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for (int l = 0; l < 3; ++l) {
						tempVector[i][j][k][l] /= outHessenbergMatrix[n - 1][n - 2];
						alertNaNOrInfinity(tempVector[i][j][k][l], "tempVector = NaN\n");
					}
				}
			}
		}
		a = scalarMultiplyLargeVectors(resultBasis[0], tempVector);
		a = scalarMultiplyLargeVectors(tempVector, tempVector);
	} else {
		printf("outHessenbergMatrix[n-1][n-2] == 0\n");
	}

	resultBasis[n - 1] = tempVector;

	return resultBasis;
}

void Simulation::generalizedMinimalResidualMethod(std::vector<MatrixElement>**** matrix, double**** rightPart, Vector3d*** outvector) {
	double norm = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart));

	if(norm == 0) {
		for(int i = 0; i < xnumber; ++i) {
			for(int j = 0; j < ynumber; ++j) {
				for(int k = 0; k < znumber; ++k) {
					outvector[i][j][k] = Vector3d(0, 0, 0);
					return;
				}
			}
		}
	}

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					rightPart[i][j][k][l] /= norm;
					for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						double value = matrix[i][j][k][l][m].value;
						matrix[i][j][k][l][m].value /= norm;
						value = matrix[i][j][k][l][m].value;
					}
				}
			}
		}
	}

	int matrixDimension = 3*(xnumber)*(ynumber)*(znumber);
	double maxError = E0.norm() / (matrixDimension*1E4);
	//double maxError = 1/(matrixDimension * 1E5);

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
	Rmatrix[1] = new double[1];
	oldRmatrix[0] = new double[1];
	oldRmatrix[1] = new double[1];

	double***** basis = new double****[1];
	basis[0] = new double***[xnumber];
	for (int i = 0; i < xnumber; ++i) {
		basis[0][i] = new double**[ynumber];
		for (int j = 0; j < ynumber; ++j) {
			basis[0][i][j] = new double*[znumber];
			for (int k = 0; k < znumber; ++k) {
				basis[0][i][j][k] = new double[3];
				for (int l = 0; l < 3; ++l) {
					basis[0][i][j][k][l] = rightPart[i][j][k][l];
				}
			}
		}
	}
	double***** newBasis;

	double**** tempResult = new double***[xnumber];
	for(int i = 0; i < xnumber; ++i) {
		tempResult[i] = new double**[ynumber];
		for(int j = 0; j < ynumber; ++j) {
			tempResult[i][j] = new double*[znumber];
			for(int k = 0; k < znumber; ++k) {
				tempResult[i][j][k] = new double[3];
			}
		}
	}

	int n = 2;
	double beta = 1.0;
	double error = beta;
	double* y = new double[1];

	double rho;
	double sigma;
	double cosn;
	double sinn;
	double module;

	double relativeError = 1;
	double maxRelativeError = 1/(matrixDimension*1E6);

	while (relativeError > maxRelativeError  && n < min2(maxGMRESIterations, matrixDimension)) {
		printf("GMRES iteration %d\n", n);
		newHessenbergMatrix = new double*[n];
		for (int i = 0; i < n; ++i) {
			newHessenbergMatrix[i] = new double[n - 1];
		}
		newBasis = arnoldiIterations(matrix, newHessenbergMatrix, n, basis, hessenbergMatrix);

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
			if(Rmatrix[i][i] > 0){
				y[i] /= Rmatrix[i][i];
			} else {
				y[i] = 0;
				printf("Rmatrix[%d][%d] = 0\n", i, i);
			}
			alertNaNOrInfinity(y[i], "y = NaN\n");
		}

		error = fabs(beta * Qmatrix[n - 1][0]);
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					outvector[i][j][k].x = 0;
					outvector[i][j][k].y = 0;
					outvector[i][j][k].z = 0;
					for (int m = 0; m < n; ++m) {
						outvector[i][j][k].x += basis[m][i][j][k][0] * y[m];
						outvector[i][j][k].y += basis[m][i][j][k][1] * y[m];
						outvector[i][j][k].z += basis[m][i][j][k][2] * y[m];
					}
				}
			}
		}

		double**** leftPart1 = multiplySpecialMatrixVector(matrix, outvector);
		double error1 = 0;
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
					for(int l = 0; l < 3; ++l){
						error1 += sqr(leftPart1[i][j][k][l] - rightPart[i][j][k][l]);
					}
					delete[] leftPart1[i][j][k];
				}
				delete[] leftPart1[i][j];
			}
			delete[] leftPart1[i];
		}
		delete[] leftPart1;
		error1 = sqrt(error1);

		double normRightPart = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart));
		relativeError = error/normRightPart;

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

	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				outvector[i][j][k].x = 0;
				outvector[i][j][k].y = 0;
				outvector[i][j][k].z = 0;
				for (int m = 0; m < n; ++m) {
					outvector[i][j][k].x += basis[m][i][j][k][0] * y[m];
					outvector[i][j][k].y += basis[m][i][j][k][1] * y[m];
					outvector[i][j][k].z += basis[m][i][j][k][2] * y[m];
					//outvector[i][j][k].x += basis[m][i][j][k][0] * y[m]*norm;
					//outvector[i][j][k].y += basis[m][i][j][k][1] * y[m]*norm;
					//outvector[i][j][k].z += basis[m][i][j][k][2] * y[m]*norm;
				}
			}
		}
	}

	double**** leftPart = multiplySpecialMatrixVector(matrix, outvector);
	error = 0;
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for(int l = 0; l < 3; ++l){
					error += sqr(leftPart[i][j][k][l] - rightPart[i][j][k][l]);
				}
				delete[] leftPart[i][j][k];
			}
			delete[] leftPart[i][j];
		}
		delete[] leftPart[i];
	}
	delete[] leftPart;

	error = sqrt(error);

	for (int i = 0; i < n; ++i) {
		delete[] Qmatrix[i];
		delete[] Rmatrix[i];
		delete[] hessenbergMatrix[i];
	}
	delete[] Qmatrix;
	delete[] Rmatrix;
	delete[] hessenbergMatrix;

	for (int m = 0; m < n; ++m) {
		for (int i = 0; i < xnumber; ++i) {
			for (int j = 0; j < ynumber; ++j) {
				for (int k = 0; k < znumber; ++k) {
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
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					result += a[i][j][k][l] * b[i][j][k][l];
				}
			}
		}
	}
	return result;
}

double Simulation::scalarMultiplyLargeVectors(Vector3d*** a, Vector3d*** b) {
	double result = 0;
	for (int i = 0; i < xnumber; ++i) {
		for (int j = 0; j < ynumber; ++j) {
			for (int k = 0; k < znumber; ++k) {
				for (int l = 0; l < 3; ++l) {
					result += a[i][j][k][l] * b[i][j][k][l];
				}
			}
		}
	}
	return result;
}