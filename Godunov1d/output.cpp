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
		double p = simulation->minP;
		double deltaP = (simulation->maxP - simulation->minP)/(pgridNumber - 1);
		for(int j = 0; j < pgridNumber; ++j){
			fprintf(distributionFile, "%17.12lf %17.12\n", p, simulation->distributionFunction[i][j]);
		}
	}
}