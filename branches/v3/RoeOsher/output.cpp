#include "output.h"
#include "constants.h"
#include "util.h"

void output(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		double t = simulation->temperatureIn(i);
		double tau = abs2(simulation->middleGrid[i]*simulation->maxRate[i]/simulation->middleVelocity[i]);
		if(simulation->middleVelocity[i] < 1){
			tau = 0;
		}
		fprintf(outFile, "%g %g %g %g %g %g %g %g %g %g", simulation->grid[i], simulation->middleVelocity[i], simulation->middleDensity[i], simulation->middlePressure[i], simulation->cosmicRayPressure[i], t, simulation->magneticInductionSum[i]-simulation->B0, simulation->magneticEnergy[i], tau, simulation->cosmicRayConcentration[i]);
		fprintf(outFile, "\n");
	}
}

void outputDistribution(FILE* distributionFile, FILE* fullDistributionFile, FILE* coordinateDistributionFile, Simulation* simulation){
	double* fullDistribution = new double[pgridNumber];
	double volume = simulation->upstreamR;
	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] = 0;
	}
	for(int i = 0; i < simulation->rgridNumber; ++i){
		int outP = injectionMomentum;
		for(int j = 0; j < pgridNumber - 1; ++j){
			double p = simulation->pgrid[j];
			fullDistribution[j] += simulation->volume(i)*simulation->distributionFunction[i][j];
		}
		fprintf(coordinateDistributionFile, "%20.10lf %g %g %g %g %g %g %g\n", simulation->grid[i], simulation->distributionFunction[i][outP], simulation->distributionFunction[i][pgridNumber-1], simulation->crflux[i][outP], simulation->crflux[i][pgridNumber-1], simulation->integratedFlux[i], simulation->distributionFunction[i][outP+5], simulation->distributionFunction[i][outP+1]);
	}

	if(simulation->shockWavePoint > 0 && simulation->shockWavePoint < simulation->rgridNumber){
		for(int j = 0; j < pgridNumber; ++j){
			fprintf(distributionFile, "%g %g %g\n", simulation->pgrid[j], simulation->distributionFunction[simulation->shockWavePoint][j], simulation->crflux[simulation->shockWavePoint][j]);
		}
	}

	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] /= volume;
		double p = simulation->pgrid[j];
		fprintf(fullDistributionFile, "%g %g\n", p, fullDistribution[j]/(cube(p)));
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
		int outP = injectionMomentum;
		for(int j = 0; j < pgridNumber - 1; ++j){
			double p = simulation->pgrid[j];
			fullDistribution[j] += simulation->volume(i)*simulation->distributionFunction[i][j];
		}
		fprintf(coordinateDistributionFile, "%20.10lf %g %g %g %g %g %g %g\n", simulation->grid[i], simulation->distributionFunction[i][outP]/(cube(simulation->pgrid[outP])), simulation->distributionFunction[i][pgridNumber-1]/(cube(simulation->pgrid[pgridNumber-1])), simulation->crflux[i][outP], simulation->crflux[i][pgridNumber-1], simulation->integratedFlux[i], simulation->distributionFunction[i][outP+5]/(cube(simulation->pgrid[outP+5])), simulation->distributionFunction[i][outP+1]/(cube(simulation->pgrid[outP+1])));
	}

	if(simulation->shockWavePoint > 0 && simulation->shockWavePoint < simulation->rgridNumber){
		for(int j = 0; j < pgridNumber; ++j){
			fprintf(distributionFile, "%g %g %g\n", simulation->pgrid[j], simulation->distributionFunction[simulation->shockWavePoint][j]/cube(simulation->pgrid[j]), simulation->crflux[simulation->shockWavePoint][j]);
		}
	}

	for(int j = 0; j < pgridNumber; ++j){
		fullDistribution[j] /= volume;
		double p = simulation->pgrid[j];
		fprintf(fullDistributionFile, "%g %g\n", p, fullDistribution[j]/(cube(p)));
	}
	delete[] fullDistribution;
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

void outputField(FILE* outFile, FILE* coordinateFile, FILE* outFull, FILE* coefFile, FILE* xfile, FILE* kfile,  Simulation* simulation){
	double* integralField = new double[kgridNumber];
	for(int k = 0; k < kgridNumber; ++k){
		integralField[k] = 0;
	}

	double volume = 0; 
	for(int i = 0; i < simulation->rgridNumber; ++i){
		fprintf(xfile, "%g\n", simulation->grid[i]);
		fprintf(coordinateFile, "%g %g\n", simulation->grid[i],simulation->magneticField[i][35]);
		if(abs(i-simulation->shockWavePoint) < 20){
			volume += simulation->volume(i);
		}
		for(int k = 0; k < kgridNumber; ++k){
			fprintf(outFull, "%g ", simulation->magneticField[i][k]);
			if(abs(i-simulation->shockWavePoint) < 20){
				integralField[k] += simulation->magneticField[i][k]*simulation->volume(i);
			}
		}
		fprintf(outFull, "\n");
	}

	for(int k = 0; k < kgridNumber; ++k){
		integralField[k] /= volume;
	}

	if(simulation->shockWavePoint > 0){
		int shockWavePoint = simulation->shockWavePoint;
		for(int k = 0; k < kgridNumber; ++k){
			fprintf(kfile, "%g\n", simulation->kgrid[k]);
			fprintf(outFile, "%g %g %g %g %g\n", simulation->kgrid[k], simulation->magneticField[shockWavePoint][k], simulation->growth_rate[shockWavePoint][k], integralField[k], abs(simulation->middleGrid[shockWavePoint]/simulation->middleVelocity[shockWavePoint])*simulation->growth_rate[shockWavePoint][k]);
		}

		for(int j = 0; j < pgridNumber; ++j){
			fprintf(coefFile, "%g %g %g %g\n", simulation->pgrid[j], simulation->diffusionCoef[shockWavePoint-1][j], simulation->diffusionCoef[shockWavePoint][j], simulation->diffusionCoef[shockWavePoint+1][j]);
		}
	}


	delete[] integralField;
}