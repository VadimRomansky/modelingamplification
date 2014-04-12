#include <time.h>
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

//изменение сетки

void Simulation::updateGrid(){
	if ((shockWavePoint < 1) || (shockWavePoint > rgridNumber - 1)) return;
	if( !shockWaveMoved) {
		return;
	}
	printf("updating grid\n");
	double shockWaveR = grid[shockWavePoint];

	double a = shockWaveR;
	double b = upstreamR - shockWaveR;
	double R1 = a;
	double R2 = b;
	double h1=0.5*rgridNumber/log(1.0+a/R1);
	double h2=(0.5*rgridNumber+1)/log(1.0+b/R2);
	for(int i=1; i < rgridNumber/2; ++ i){ 
		tempGrid[i] = shockWaveR + R1*(1 - exp(-(1.0*(i+1)-0.5*rgridNumber)/h1));
	}
	for(int i=rgridNumber/2; i < rgridNumber; ++i){
		tempGrid[i] = shockWaveR + R2*(exp((1.0*(i+1)-0.5*rgridNumber)/h2)-1.0);
	}
	tempGrid[rgridNumber] = upstreamR;
	for(int i = 1; i <= rgridNumber; ++i){
		if(tempGrid[i] <= tempGrid[i-1]){
			printf("tempGrid[i] <= tempGrid[i-1]\n");
		}
	}

	redistributeValues();
}


//перераспределение величин между ячейками новой сетки

void Simulation::redistributeValues(){
	int oldCount = 1;
	double tempDensity = 0;
	double tempMomentum = 0;
	double tempEnergy = 0;
	double* newDensity = new double[rgridNumber];
	double* newMomentum = new double[rgridNumber];
	double* newEnergy = new double[rgridNumber];
	double** newDistributionFunction = new double*[rgridNumber];
	double* tempFunction = new double[pgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		newDistributionFunction[i] = new double[pgridNumber];
		double newVolume = 4*pi*(cube(tempGrid[i+1]) - cube(tempGrid[i]))/3;
		bool oldChanged = false;
		if(tempGrid[i+1] > grid[oldCount]){
			tempDensity = middleDensity[oldCount - 1]*4*pi*(cube(grid[oldCount]) - cube(tempGrid[i]))/(3*newVolume);
			tempMomentum = momentum(oldCount - 1)*4*pi*(cube(grid[oldCount]) - cube(tempGrid[i]))/(3*newVolume);
			tempEnergy = energy(oldCount - 1)*4*pi*(cube(grid[oldCount]) - cube(tempGrid[i]))/(3*newVolume);
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j]*4*pi*(cube(grid[oldCount]) - cube(tempGrid[i]))/(3*newVolume);
			}
			++oldCount;
			oldChanged = true;
		} else {
			tempDensity = middleDensity[oldCount - 1];
			tempMomentum = momentum(oldCount - 1);
			tempEnergy = energy(oldCount - 1);
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j];
			}
		}
		while(grid[oldCount] < tempGrid[i+1]){
			tempDensity += middleDensity[oldCount - 1]*volume(oldCount - 1)/newVolume;
			tempMomentum += momentum(oldCount - 1)*volume(oldCount - 1)/newVolume;
			tempEnergy += energy(oldCount - 1)*volume(oldCount - 1)/newVolume;
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j]*volume(oldCount - 1)/newVolume;
			}
			++oldCount;
			oldChanged = true;
			if(oldCount > rgridNumber + 1){
				printf("oldCount > rgridNUmber + 1\n");
			}
		}
		if(oldChanged){
			tempDensity += middleDensity[oldCount - 1]*4*pi*(cube(tempGrid[i+1]) - cube(grid[oldCount - 1]))/(3*newVolume);
			tempMomentum += momentum(oldCount - 1)*4*pi*(cube(tempGrid[i+1]) - cube(grid[oldCount - 1]))/(3*newVolume);
			tempEnergy += energy(oldCount - 1)*4*pi*(cube(tempGrid[i+1]) - cube(grid[oldCount - 1]))/(3*newVolume);
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j]*4*pi*(cube(tempGrid[i+1]) - cube(grid[oldCount - 1]))/(3*newVolume);
			}
		}
		newDensity[i] = tempDensity;
		alertNaNOrInfinity(newDensity[i], "newDensity = NaN");
		alertNegative(newDensity[i], "newDensity < 0");
		newMomentum[i] = tempMomentum;
		alertNaNOrInfinity(newMomentum[i], "newMomentum = NaN");
		newEnergy[i] = tempEnergy;
		alertNaNOrInfinity(newEnergy[i], "newEnergy = NaN");
		alertNegative(newEnergy[i], "newEnergy < 0");
		for(int j = 0; j < pgridNumber; ++j){
			newDistributionFunction[i][j] = tempFunction[j];
			alertNaNOrInfinity(newDistributionFunction[i][j], "newDistribution = NaN");
			alertNegative(newDistributionFunction[i][j], "newDistribution < 0");
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		grid[i + 1] = tempGrid[i + 1];
		gridsquare[i+1] = grid[i+1]*grid[i+1];
		middleGrid[i] = (grid[i] + grid[i+1])/2;
		volumeFactor[i] = (cube(grid[i+1]) - cube(grid[i]))/3;
		deltaR[i] = grid[i+1] - grid[i];
		if(i == 0){
			middleDeltaR[i] = middleGrid[i];
		} else {
			middleDeltaR[i] = middleGrid[i] - middleGrid[i-1];
		}

		middleDensity[i] = newDensity[i];
		if(newDensity[i] <= epsilon*density0){
			middleVelocity[i] = 0;
		} else {
			middleVelocity[i] = newMomentum[i]/newDensity[i];
		}
		double tempPressure = (newEnergy[i] - middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2)*(gamma - 1);
		
		if(tempPressure < 0){
			middlePressure[i] = 0.01*min2(middlePressure[i+1],middlePressure[i]);
			printf("pressure < 0\n");
		} else {
			middlePressure[i] = tempPressure;
		}
		alertNegative(middlePressure[i], "middlePressure < 0");
		alertNaNOrInfinity(middlePressure[i], "middlePressure = NaN");

		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = newDistributionFunction[i][j];
		}
	}
	if(newDensity[rgridNumber - 1] < middleDensity[rgridNumber - 1]){
		printf("aaa\n");
	}

	for(int i = 0; i < rgridNumber; ++i){
		delete[] newDistributionFunction[i];
	}
	delete[] newDistributionFunction;
	delete[] newDensity;
	delete[] newMomentum;
	delete[] newEnergy;
	delete[] tempFunction;
}