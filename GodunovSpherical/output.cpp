#include "stdafx.h"
#include "output.h"
#include "constants.h"

void output(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		double t = simulation->temperatureIn(i);
		//fprintf(outFile,"%lf %lf %lf %lf %lf\n", simulation->grid[i], simulation->middleVelocity[i], simulation->middleDensity[i], simulation->middlePressure[i], t);
	    fprintf(outFile,"%17.12lf %17.12lf %38.30lf %28.20lf %17.12lf\n", simulation->grid[i], simulation->middleVelocity[i], simulation->middleDensity[i], simulation->middlePressure[i], t);
		//fprintf(outFile,"%lf %lf %lf %lf %lf\n", simulation->grid[i], simulation->middleVelocity[i], simulation->middleDensity[i], simulation->middlePressure[i], simulation->temperatureIn(i));
	}
}

void outputDistribution(FILE* distributionFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		for(int j = 0; j < pgridNumber; ++j){
			double p = (simulation->pgrid[j+1] + simulation->pgrid[j])/2;
			fprintf(distributionFile, "%20.15lf %20.15lf\n", p, simulation->distributionFunction[i][j]);
		}
	}
}