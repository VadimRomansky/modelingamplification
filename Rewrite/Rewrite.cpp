// Rewrite.cpp : Defines the entry point for the console application.
//

#include "stdio.h"

const double massProton = 1.67262177E-24;
const double kBoltzman = 1.3806488E-16;
const double speed_of_light = 2.99792458E10;
const double electron_charge = 0.0000000004803;
const double pi = 3.1415926535897932384626433832795028841971693993751;


int main()
{
	int rgridNumber = 1000;
	int pgridNumber = 150;
	int kgridNumber = 100;
	int iterations = 11;
	int count = 11;

	FILE* profile = fopen("./input/tamc_radial_profile.dat","r");
	double** all = new double*[iterations*rgridNumber];
	for(int i = 0; i < iterations*rgridNumber; ++i){
		all[i] = new double[count];
	}

	for(int i = 0; i < iterations*rgridNumber; ++i){
		for(int j = 0; j < count; ++j){
			fscanf(profile,"%lf",&all[i][j]);
		}
	}


	double** velocity = new double*[rgridNumber];
	double** density = new double*[rgridNumber];
	double** pressure = new double*[rgridNumber];
	double** temperature = new double*[rgridNumber];
	double** field = new double*[rgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		velocity[i] = new double[4];
		density[i] = new double[4];
		pressure[i] = new double[4];
		temperature[i] = new double[4];
		field[i] = new double[4];
	}

	int a = iterations - 1;
	int b = 3*a/4;
	int c = a/2;

	for(int i = 0; i < rgridNumber; ++i){
		velocity[i][0] = all[a*rgridNumber + i][0];
		velocity[i][1] = all[a*rgridNumber + i][1];
		velocity[i][2] = all[b*rgridNumber + i][1];
		velocity[i][3] = all[c*rgridNumber + i][1];

		density[i][0] = all[a*rgridNumber + i][0];
		density[i][1] = all[a*rgridNumber + i][2];
		density[i][2] = all[b*rgridNumber + i][2];
		density[i][3] = all[c*rgridNumber + i][2];

		pressure[i][0] = all[a*rgridNumber + i][0];
		pressure[i][1] = all[a*rgridNumber + i][3];
		pressure[i][2] = all[b*rgridNumber + i][3];
		pressure[i][3] = all[c*rgridNumber + i][3];

		temperature[i][0] = all[a*rgridNumber + i][0];
		temperature[i][1] = all[a*rgridNumber + i][5];
		temperature[i][2] = all[b*rgridNumber + i][5];
		temperature[i][3] = all[c*rgridNumber + i][5];

		field[i][0] = all[a*rgridNumber + i][0];
		field[i][1] = all[a*rgridNumber + i][6];
		field[i][2] = all[b*rgridNumber + i][6];
		field[i][3] = all[c*rgridNumber + i][6];
	}
	fclose(profile);

	FILE* vel = fopen("./output/velocity.dat","w");
	FILE* rho = fopen("./output/density.dat","w");
	FILE* p = fopen("./output/pressure.dat","w");
	FILE* t = fopen("./output/temperature.dat","w");
	FILE* B = fopen("./output/field.dat","w");
	for(int i = 0; i < rgridNumber; ++i){
		fprintf(vel, "%g %g %g %g\n",velocity[i][0], velocity[i][1], velocity[i][2], velocity[i][3]);
		fprintf(rho, "%g %g %g %g\n",density[i][0], density[i][1], density[i][2], density[i][3]);
		fprintf(p, "%g %g %g %g\n",pressure[i][0], pressure[i][1], pressure[i][2], pressure[i][3]);
		fprintf(t, "%g %g %g %g\n",temperature[i][0], temperature[i][1], temperature[i][2], temperature[i][3]);
		fprintf(B, "%g %g %g %g\n",field[i][0], field[i][1], field[i][2], field[i][3]);
	}
	fclose(vel);
	fclose(rho);
	fclose(p);
	fclose(t);
	fclose(B);


	for(int i = 0; i < iterations*rgridNumber; ++i){
		delete[] all[i];
	}
	delete[] all;

	for(int i = 0; i < rgridNumber; ++i){
		delete[] velocity[i];
		delete[] density[i];
		delete[] pressure[i];
		delete[] temperature[i];
		delete[] field[i];
	}

	FILE* distributionIn = fopen("./input/distribution.dat","r");
	double** distributionAll = new double*[pgridNumber*iterations];
	for(int i = 0; i < pgridNumber*iterations; ++i){
		distributionAll[i] = new double[3];
	}
	double** distribution = new double*[pgridNumber];
	for(int i = 0; i < pgridNumber; ++i){
		distribution[i] = new double[4];
	}

	for(int i = 0; i < pgridNumber*iterations; ++i){
		for(int j = 0; j < 3; ++j){
			fscanf(distributionIn, "%lf", &distributionAll[i][j]);
		}
	}
	fclose(distributionIn);

	for(int i = 0; i < pgridNumber; ++i){
		distribution[i][0] = distributionAll[a*pgridNumber + i][0]/(massProton*speed_of_light);
		double p4 = distribution[i][0]*distribution[i][0]*distribution[i][0]*distribution[i][0];
		distribution[i][1] = distributionAll[a*pgridNumber + i][1]*p4;
		distribution[i][2] = distributionAll[b*pgridNumber + i][1]*p4;
		distribution[i][3] = distributionAll[c*pgridNumber + i][1]*p4;
	}
	for(int i = 0; i < pgridNumber*iterations; ++i){
		delete[] distributionAll[i];
	}
	delete[] distributionAll;

	FILE* distributionOut = fopen("./output/distribution.dat","w");

	for(int i = 0; i < pgridNumber; ++i){
		fprintf(distributionOut, "%g %g %g %g\n", distribution[i][0], distribution[i][1], distribution[i][2], distribution[i][3]);
	}

	fclose(distributionOut);

	//for(int i = 0; i < pgridNumber; ++i){
		//delete[] distribution[i];
	//}
	//delete[] distribution;

	return 0;
}

