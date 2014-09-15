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
				matrix[i][j] = 10;
			} else {
				matrix[i][j] = 1;
			}
		}
	}
	matrix[1][2] = 2;
	matrix[5][3] = 2;
	matrix[2][4] = 2;

	for(int i = 0; i < number; ++i){
		rightPart[i] = 0;
	}
	rightPart[0] = 1/sqrt(2.0);
	rightPart[1] = 1/sqrt(2.0);

	generalizedMinimalResidualMethod(matrix, rightPart);


	return 0;
}

