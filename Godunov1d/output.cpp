#include "output.h"
#include "constants.h"
#include "util.h"

void output(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		double t = simulation->temperatureIn(i);
		//fprintf(outFile,"%lf %lf %lf %lf %lf\n", simulation->grid[i], simulation->middleVelocity[i], simulation->middleDensity[i], simulation->middlePressure[i], t);
	    fprintf(outFile,"%17.12lf %17.12lf %38.30lf %28.20lf %28.20lf %17.12lf\n", simulation->grid[i], simulation->middleVelocity[i], simulation->middleDensity[i], simulation->middlePressure[i], simulation->cosmicRayPressure[i], t);
		//fprintf(outFile,"%lf %lf %lf %lf %lf\n", simulation->grid[i], simulation->middleVelocity[i], simulation->middleDensity[i], simulation->middlePressure[i], simulation->temperatureIn(i));
	}
}

void outputDistribution(FILE* distributionFile, FILE* fullDistributionFile, FILE* coordinateDistributionFile, Simulation* simulation){
	double* fullDistribution = new double[pgridNumber];
	double volume = simulation->upstreamR;
	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] = 0;
	}
	for(int i = 0; i < simulation->rgridNumber; ++i){
		for(int j = 0; j < pgridNumber - 1; ++j){
			//double p = (simulation->pgrid[j+1] + simulation->pgrid[j])/2;
			double p = simulation->pgrid[j];
			//fprintf(distributionFile, "%30.20lf %g\n", p, simulation->distributionFunction[i][j]);
			fullDistribution[j] += simulation->volume(i)*simulation->distributionFunction[i][j];
		}
		//fprintf(coordinateDistributionFile, "%20.10lf %g\n", simulation->grid[i], simulation->distributionFunction[i][injectionMomentum]/(cube(simulation->pgrid[injectionMomentum])));
		fprintf(coordinateDistributionFile, "%20.10lf %g\n", simulation->grid[i], simulation->distributionFunction[i][injectionMomentum]);
	}

	if(simulation->shockWavePoint > 0 && simulation->shockWavePoint < simulation->rgridNumber){
		for(int j = 0; j < pgridNumber; ++j){
			//fprintf(distributionFile, "%g %g\n", simulation->pgrid[j], simulation->distributionFunction[simulation->shockWavePoint][j]/cube(simulation->pgrid[j]));
			fprintf(distributionFile, "%g %g\n", simulation->pgrid[j], simulation->distributionFunction[simulation->shockWavePoint][j]);
		}
	}

	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] /= volume;
		double p = simulation->pgrid[j];
		//fprintf(fullDistributionFile, "%g %g\n", p, fullDistribution[j]/(cube(p)));
		fprintf(fullDistributionFile, "%g %g\n", p, fullDistribution[j]);
	}
	delete[] fullDistribution;
}

void outputNewGrid(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		fprintf(outFile, "%d %17.12lf\n", i, simulation->tempGrid[i]);
	}
}