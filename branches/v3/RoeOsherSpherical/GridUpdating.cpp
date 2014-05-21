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
	double rightR = grid[rgridNumber] - shockWaveR;

	int leftPoints = shockWavePoint- 4;
	int rightPonts = rgridNumber - 3 - shockWavePoint + 1;

	double R1 = grid[shockWavePoint-4];
	double R2 = grid[rgridNumber] - grid[shockWavePoint+3];
	double a = 10000;
	double b = 10000;
	double h1 = leftPoints/log(1.0+a);
	double h2 = rightPonts/log(1.0+b);

	for(int i = 0; i < shockWavePoint - 3; ++i){
		tempGrid[i] = (R1/a)*(1 - exp(-(1.0*(i+1)-leftPoints)/h1)) + R1;
	}

	for(int i = shockWavePoint -3; i < shockWavePoint + 3; ++i){
		tempGrid[i] = grid[i];
	}

	for(int i = shockWavePoint + 3; i < rgridNumber; ++i){
		tempGrid[i] = (R2/b)*(exp((1.0*(i+1)-rgridNumber)/h2)-1.0) + upstreamR;
	}

	for(int i = 1; i <= rgridNumber; ++i){
		if(tempGrid[i] < tempGrid[i-1]){
			printf("grid[i] < grid[i-1]\n");
		}
	}

	tempGrid[0] = 0;
	tempGrid[rgridNumber] = grid[rgridNumber];

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
		double newVolume = tempGrid[i+1] - tempGrid[i];
		bool oldChanged = false;
		if(tempGrid[i+1] > grid[oldCount]){
			tempDensity = middleDensity[oldCount - 1]*(grid[oldCount] - tempGrid[i])/(newVolume);
			tempMomentum = momentum(oldCount - 1)*(grid[oldCount] - tempGrid[i])/(newVolume);
			tempEnergy = energy(oldCount - 1)*(grid[oldCount] - tempGrid[i])/(newVolume);
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j]*(grid[oldCount] - tempGrid[i])/(newVolume);
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
			tempDensity += middleDensity[oldCount - 1]*(tempGrid[i+1] - grid[oldCount - 1])/(newVolume);
			tempMomentum += momentum(oldCount - 1)*(tempGrid[i+1] - grid[oldCount - 1])/(newVolume);
			tempEnergy += energy(oldCount - 1)*(tempGrid[i+1] - grid[oldCount - 1])/(newVolume);
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j]*(tempGrid[i+1] - grid[oldCount - 1])/(newVolume);
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
		middleGrid[i] = (grid[i] + grid[i+1])/2;
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