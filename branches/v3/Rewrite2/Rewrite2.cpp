// Rewrite2.cpp : Defines the entry point for the console application.
//

#include "stdio.h"

const double massProton = 1.67262177E-24;
const double kBoltzman = 1.3806488E-16;
const double speed_of_light = 2.99792458E10;
const double electron_charge = 0.0000000004803;
const double pi = 3.1415926535897932384626433832795028841971693993751;

int main()
{
	int pgridNumber = 150;
	int kgridNumber = 100;
	int iterations = 72;

	FILE* distributionIn1 = fopen("./input1/distribution.dat","r");
	FILE* distributionIn2 = fopen("./input2/distribution.dat","r");
	double** distributionAll1 = new double*[pgridNumber*iterations];
	double** distributionAll2 = new double*[pgridNumber*iterations];
	for(int i = 0; i < pgridNumber*iterations; ++i){
		distributionAll1[i] = new double[3];
		distributionAll2[i] = new double[3];
	}
	double** distribution = new double*[pgridNumber];
	for(int i = 0; i < pgridNumber; ++i){
		distribution[i] = new double[3];
	}

	for(int i = 0; i < pgridNumber*iterations; ++i){
		for(int j = 0; j < 3; ++j){
			fscanf(distributionIn1, "%lf", &distributionAll1[i][j]);
			fscanf(distributionIn2, "%lf", &distributionAll2[i][j]);
		}
	}
	fclose(distributionIn1);
	fclose(distributionIn2);

	int a = iterations - 1;

	for(int i = 0; i < pgridNumber; ++i){
		distribution[i][0] = distributionAll1[a*pgridNumber + i][0]/(massProton*speed_of_light);
		double p4 = distribution[i][0]*distribution[i][0]*distribution[i][0]*distribution[i][0];
		distribution[i][1] = distributionAll1[(a-0)*pgridNumber + i][1]*p4;
		distribution[i][2] = distributionAll2[a*pgridNumber + i][1]*p4;
	}
	for(int i = 0; i < pgridNumber*iterations; ++i){
		delete[] distributionAll1[i];
		delete[] distributionAll2[i];
	}
	delete[] distributionAll1;
	delete[] distributionAll2;

	FILE* distributionOut = fopen("./output/distribution2.dat","w");

	for(int i = 0; i < pgridNumber; ++i){
		fprintf(distributionOut, "%g %g %g\n", distribution[i][0], distribution[i][1], distribution[i][2]);
	}

	fclose(distributionOut);

	return 0;
}

