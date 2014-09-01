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
				matrix[i][j] = 0;
			}
		}
	}

	matrix[0][3] = 1;
	matrix[1][5] = 1;
	matrix[6][7] = 1;

	for(int i = 0; i < number; ++i){
		rightPart[i] = 1;
	}

	generalizedMinimalResidualMethod(matrix, rightPart);


	return 0;
}

