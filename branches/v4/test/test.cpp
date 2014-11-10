// test.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "specialmath.h"


int main()
{
	double** matrix = new double*[number];
	double* rightPart = new double[number];

	for(int i = 0; i < number; ++i){
		matrix[i] = new double[number];
		for(int j = 0; j < number; ++j){
			if(i == j){
				matrix[i][j] = 0.01;
			} else {
				matrix[i][j] = 0;
			}
		}
	}
	matrix[0][1] = 1;
	matrix[1][2] = 1;
	matrix[2][3] = 1;
	matrix[3][4] = 1;
	matrix[4][5] = 1;
	matrix[5][6] = 1;

	for(int i = 0; i < number; ++i){
		rightPart[i] = 0;
	}
	rightPart[number-1] = 1;

	/*for(int i = 0; i < number; ++i) {
		matrix[i] = new double[number];
		for(int j = 0; j < number; ++j) {
			matrix[i][j] = 0;
		}
		rightPart[i] = 1;
	}*/

	for(int i = 0; i < number; ++i) {
		for(int j = 0; j < number; ++j) {
			printf("%lf     ", matrix[i][j]);
		}
		printf("   %lf\n", rightPart[i]);
	}

	double* result = generalizedMinimalResidualMethod(matrix, rightPart);

	printf("result\n");
	for(int i = 0; i < number; ++i) {
		printf("%15.10lf\n", result[i]);
	}


	return 0;
}
