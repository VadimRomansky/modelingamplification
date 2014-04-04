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
			double p = simulation->pgrid[j];
			fullDistribution[j] += simulation->volume(i)*simulation->distributionFunction[i][j];
		}
		fprintf(coordinateDistributionFile, "%20.10lf %g\n", simulation->grid[i], simulation->distributionFunction[i][injectionMomentum]);
	}

	if(simulation->shockWavePoint > 0 && simulation->shockWavePoint < simulation->rgridNumber){
		for(int j = 0; j < pgridNumber; ++j){
			fprintf(distributionFile, "%g %g\n", simulation->pgrid[j], simulation->distributionFunction[simulation->shockWavePoint][j]);
		}
	}

	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] /= volume;
		double p = simulation->pgrid[j];
		fprintf(fullDistributionFile, "%g %g\n", p, fullDistribution[j]);
	}
	delete[] fullDistribution;
}

void outputDistributionP3(FILE* distributionFile, FILE* fullDistributionFile, FILE* coordinateDistributionFile, Simulation* simulation){
	double* fullDistribution = new double[pgridNumber];
	double volume = simulation->upstreamR;
	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] = 0;
	}
	for(int i = 0; i < simulation->rgridNumber; ++i){
		for(int j = 0; j < pgridNumber - 1; ++j){
			double p = simulation->pgrid[j];
			fullDistribution[j] += simulation->volume(i)*simulation->distributionFunction[i][j];
		}
		fprintf(coordinateDistributionFile, "%20.10lf %g\n", simulation->grid[i], simulation->distributionFunction[i][injectionMomentum]/(cube(simulation->pgrid[injectionMomentum])));
	}

	if(simulation->shockWavePoint > 0 && simulation->shockWavePoint < simulation->rgridNumber){
		for(int j = 0; j < pgridNumber; ++j){
			fprintf(distributionFile, "%g %g\n", simulation->pgrid[j], simulation->distributionFunction[simulation->shockWavePoint][j]/cube(simulation->pgrid[j]));
		}
	}

	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] /= volume;
		double p = simulation->pgrid[j];
		fprintf(fullDistributionFile, "%g %g\n", p, fullDistribution[j]/(cube(p)));
	}
	delete[] fullDistribution;
}

void outputDerivativeForDebug(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < pgridNumber; ++i){
		fprintf(outFile, "%g %g %g\n", simulation->pgrid[i], simulation->distrFunDerivative[i], simulation->distrFunDerivative2[i]);
	}
}

void outputNewGrid(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		fprintf(outFile, "%d %17.12lf\n", i, simulation->tempGrid[i]);
	}
}

void outMatrix(double* a, double* c, double* b, int N, double* f, double* x){
	FILE* file = fopen("output/matrix.dat","w");
	fprintf(file, "%g %g %g %g %g\n" , 0.0, c[0], 0.0, f[0], x[0]);
	for(int i = 1; i <= N; ++i){
		fprintf(file, "%g %g %g %g %g\n" , a[i-1], c[i], b[i-1], f[i], x[i]);
	}
	fclose(file);
}