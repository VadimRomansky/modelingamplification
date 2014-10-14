#ifndef _SPECIAL_MATH_H_
#define _SPECIAL_MATH_H_

const int number = 3;

double** arnoldiIterations(double** matrix, double** outHessenbergMatrix, int n, double** prevBasis, double** prevHessenbergMatrix);
void generalizedMinimalResidualMethod(double** matrix, double* rightPart);
double scalarMultiplyLargeVectors(double* a, double* b);
double* multiplyMatrixVector(double** matrix, double* vector);

#endif