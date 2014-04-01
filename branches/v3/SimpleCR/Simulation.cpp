#include <list>
#include <time.h>
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

//конструктор
Simulation::Simulation(){
	myTime = 0;
	shockWavePoint = -1;
	injectedParticles = 0;
}

//деструктор
Simulation::~Simulation(){
	delete[] pgrid;
	delete[] logPgrid;
	delete[] grid;
	delete[] middleGrid;
	delete[] deltaR;
	delete[] middleDeltaR;


	for(int i = 0; i < rgridNumber; ++i){
		delete[] distributionFunction[i];
		delete[] tempDistributionFunction[i];
	}

	delete[] distributionFunction;
	delete[] tempDistributionFunction;
}

void Simulation::initializeProfile(){
	downstreamR = -upstreamR;

	pgrid = new double[pgridNumber];
	logPgrid = new double[pgridNumber];
	grid = new double[rgridNumber + 1];
	middleGrid = new double[rgridNumber];
	deltaR = new double[rgridNumber];
	middleDeltaR = new double[rgridNumber];

	distributionFunction = new double*[rgridNumber + 1];
	tempDistributionFunction = new double*[rgridNumber + 1];
	for(int i = 0; i < rgridNumber + 1; i ++){
		distributionFunction[i] = new double[pgridNumber];
		tempDistributionFunction[i] = new double[pgridNumber];
		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = 0;
		}
	}

	double r = downstreamR;
	double pressure0 = density0*kBoltzman*temperature/massProton;

	minP = massProton*speed_of_light/10;
	maxP = minP*100000000;

	deltaR0 = (upstreamR - downstreamR)/rgridNumber;
	for(int i = 0; i < rgridNumber + 1; ++i){
		grid[i] = r;
		middleGrid[i] = r + deltaR0/2;
		deltaR[i] = deltaR0;
		r += deltaR0;
	}

	logPgrid[0] = log(minP);
	logPgrid[pgridNumber - 1] = log(maxP);
	deltaLogP = (logPgrid[pgridNumber - 1] - logPgrid[0])/(pgridNumber - 1);
	pgrid[0] = minP;
	for(int i = 1; i < pgridNumber; ++i){
		logPgrid[i] = logPgrid[i-1] + deltaLogP;
		pgrid[i] = exp(logPgrid[i]);
	}
	pgrid[pgridNumber-1] = maxP;

	for(int i = 0; i < rgridNumber; ++i){
		if(i == 0){
			middleDeltaR[i] = middleGrid[i];
		} else {
			middleDeltaR[i] = middleGrid[i] - middleGrid[i-1];
		}
		double p;
		for(int j = 0; j < pgridNumber; ++j){
			p = exp(logPgrid[j]);
			distributionFunction[i][j] = 0.1*density(grid[i])*minP*minP/(massProton*maxP*pgrid[j]*pgrid[j]*pgrid[j]*pgrid[j]);
			//double x = -(sqrt(sqr(massProton*speed_of_light*speed_of_light) + sqr(p*speed_of_light))-massProton*speed_of_light*speed_of_light)/(kBoltzman*temperatureIn(i));
			//distributionFunction[i][j] = exp(x);
			//distributionFunction[i][j] = epsilon;
		}
	}

	grid[rgridNumber] = upstreamR;
}

//главная функция
void Simulation::simulate(){
	printf("initialization\n");

	initializeProfile();
	updateParameters();

	printf("creating files\n");
	FILE* outDistribution;
	FILE* outFullDistribution;
	FILE* outCoordinateDistribution;
	fopen_s(&outDistribution, "./output/distribution.dat","w");
	fopen_s(&outFullDistribution, "./output/fullDistribution.dat","w");
	fopen_s(&outCoordinateDistribution, "./output/coordinateDistribution.dat","w");
	outputDistribution(outDistribution, outFullDistribution, outCoordinateDistribution, this);
	fclose(outCoordinateDistribution);
	fclose(outFullDistribution);
	fclose(outDistribution);

	deltaT = 5000;

	currentIteration = 0;
	//основной цикл
	while(myTime < maxTime && currentIteration < iterationNumber){
		++currentIteration;
		printf("iteration № %d\n", currentIteration);
		printf("time = %lf\n", myTime);
		printf("solving\n");

		evaluateCR();


		myTime = myTime + deltaT;

		updateParameters();
		if(currentIteration % writeParameter == 0){
			//вывод на некоторых итерациях

			fopen_s(&outDistribution, "./output/distribution.dat","a");
			fopen_s(&outFullDistribution, "./output/fullDistribution.dat","a");
			fopen_s(&outCoordinateDistribution, "./output/coordinateDistribution.dat","a");

			outputDistributionP3(outDistribution, outFullDistribution, outCoordinateDistribution, this);

			fclose(outCoordinateDistribution);
			fclose(outFullDistribution);
			fclose(outDistribution);

		}
	}
}

double Simulation::velocity(double x){
	if(x <= 0){
		return U0;
	} else {
		return U0/4;
	}
}

double Simulation::density(double x){
	if(x <= 0){
		return density0/4;
	} else {
		return density0;
	}
}

double Simulation::volume(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i <= rgridNumber) {
		return deltaR[i];
	} else {
		printf("i > rgridNumber");
	}
	return 0;
}

//подсчет полной массы энергии и импульса

void Simulation::updateParameters(){
	totalParticles = 0;
	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < pgridNumber; ++j){
			double dp;
			if(j == 0){
				dp = (pgrid[j + 1] - pgrid[j]);
			} else if(j == pgridNumber -1){
				dp = (pgrid[j] - pgrid[j - 1]);
			} else {
				dp = (pgrid[j + 1] - pgrid[j - 1])/2;
			}
			totalParticles += distributionFunction[i][j]*volume(i)*dp/pgrid[j];
		}
	}
}
