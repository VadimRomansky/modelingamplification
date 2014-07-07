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
	double tempRightGridLevel = (0.5*rightR - minDeltaR)/(0.5*rightR - deltaR0);
	if(tempRightGridLevel < 1){
		return;
	}

	int rightPointsExp = log(deltaR0/minDeltaR)/log(tempRightGridLevel);
	if(rightPointsExp > 6*rgridNumber/10){
		rightPointsExp = 6*rgridNumber/10;

	/*int rightPointsExp = log(deltaR0/minDeltaR)/log(tempRightGridLevel);
	if(rightPointsExp > 7*rgridNumber/10){
		rightPointsExp = 7*rgridNumber/10;*/

	}
	if(rightPointsExp <= 1){
		printf("rightPoints <= 1!!!\n");
	}
	if(rightPointsExp < rgridNumber/10){
		return;
	}
	int rightPointsLinear = rightPointsExp/5;

	//double deltaR0Left = 0.5*shockWaveR/(rgridNumber - rightPointsExp - rightPointsLinear);
	//double tempLeftGridLevel = (0.5*shockWaveR - minDeltaR)/(0.5*shockWaveR - deltaR0Left);
	//int leftPointsExp = log(deltaR0Left/minDeltaR)/log(tempLeftGridLevel);
	//leftPointsExp = min2((rgridNumber - rightPointsExp - rightPointsLinear)/10, leftPointsExp);
	// leftPointsLinear = rgridNumber - leftPointsExp - rightPointsExp - rightPointsLinear;
	int leftPointsLinear = rgridNumber - rightPointsExp - rightPointsLinear;

	tempGrid[0] = 0;
	tempGrid[rgridNumber] = grid[rgridNumber];
	//tempGrid[leftPointsLinear + leftPointsExp] = shockWaveR;
	tempGrid[leftPointsLinear] = shockWaveR;
	//shockWavePoint = leftPointsLinear + leftPointsExp;
	shockWavePoint = leftPointsLinear;

	/*double dRLeft = deltaR0Left*exp(-leftPointsExp*log(tempLeftGridLevel));
	if(dRLeft < deltaR0Left){
		printf("dRLeft < minDeltaR");
	}
	double logLevelLeft = log(tempLeftGridLevel);
	for(int i =  leftPointsLinear + leftPointsExp - 1; i > leftPointsLinear; --i){
		tempGrid[i] = tempGrid[i + 1] - dRLeft*exp(-(i - leftPointsLinear - leftPointsExp + 1)*logLevelLeft);
	}*/

	//double leftDeltaR = (tempGrid[leftPointsLinear + 1] - grid[0])/(leftPointsLinear + 1);

	double leftDeltaR = (shockWaveR - grid[0])/(leftPointsLinear + 1);

	for(int i = 1; i < leftPointsLinear; ++i){
		tempGrid[i] = tempGrid[i-1] + leftDeltaR;
	}
	double dR = deltaR0*exp(-rightPointsExp*log(tempRightGridLevel));
	double logLevel = log(tempRightGridLevel);
	/*for(int i = leftPointsLinear + leftPointsExp + 1; i < leftPointsLinear + leftPointsExp + rightPointsExp; ++i){
		tempGrid[i] = tempGrid[i - 1] + dR*exp((i - leftPointsLinear - leftPointsExp - 1)*logLevel);
	}*/
	//dR = (grid[rgridNumber] - tempGrid[leftPointsLinear + rightPointsExp - 1])/(rightPointsLinear + 1);
	for(int i = leftPointsLinear + 1; i < leftPointsLinear + rightPointsExp; ++i){
		tempGrid[i] = tempGrid[i - 1] + dR*exp((i - leftPointsLinear - 1)*logLevel);
	}
	dR = (grid[rgridNumber] - tempGrid[leftPointsLinear + rightPointsExp - 1])/(rightPointsLinear + 1);
	for(int i = leftPointsLinear + rightPointsExp; i < rgridNumber; ++i){
		tempGrid[i] = tempGrid[i - 1] + dR;
		if(tempGrid[i] > grid[rgridNumber]){
			printf("grid > upstreamR\n");
		}
	}
	tempGrid[rgridNumber - 1] = (tempGrid[rgridNumber - 2] + tempGrid[rgridNumber])/2;
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
		double tempPressure = (newEnergy[i] - middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2)*(_gamma - 1);
		
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