#include "stdafx.h"
#include "output.h"
#include "constants.h"
#include "util.h"

void output(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		double t = simulation->temperatureIn(i);
		//fprintf(outFile,"%lf %lf %lf %lf %lf\n", simulation->grid[i], simulation->middleVelocity[i], simulation->middleDensity[i], simulation->middlePressure[i], t);
	    fprintf(outFile,"%17.12lf %17.12lf %38.30lf %28.20lf %17.12lf\n", simulation->grid[i], simulation->middleVelocity[i], simulation->middleDensity[i], simulation->middlePressure[i], t);
		//fprintf(outFile,"%lf %lf %lf %lf %lf\n", simulation->grid[i], simulation->middleVelocity[i], simulation->middleDensity[i], simulation->middlePressure[i], simulation->temperatureIn(i));
	}
}

void outputDistribution(FILE* distributionFile, FILE* fullDistributionFile, Simulation* simulation){
	double* fullDistribution = new double[pgridNumber];
	double volume = 4*pi*cube(simulation->upstreamR);
	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] = 0;
	}
	for(int i = 0; i < simulation->rgridNumber; ++i){
		for(int j = 0; j < pgridNumber - 1; ++j){
			double p = (simulation->pgrid[j+1] + simulation->pgrid[j])/2;
			fprintf(distributionFile, "%30.20lf %30.20lf\n", p, simulation->distributionFunction[i][j]);
			fullDistribution[j] += simulation->volume(i)*simulation->distributionFunction[i][j];
		}
	}

	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] /= volume;
		double p = (simulation->pgrid[j+1] + simulation->pgrid[j])/2;
		fprintf(fullDistributionFile, "%30.20lf %30.20lf\n", p, fullDistribution[j]);
	}
	delete[] fullDistribution;
}