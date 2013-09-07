#include "stdafx.h"
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

Simulation::Simulation(){
	cycleBound = true;
	time = 0;
}

Simulation::~Simulation(){
	delete[] grid;
	delete[] pointDensityLeft;
	delete[] pointDensityRight;
	delete[] pointVelocity;
	delete[] pointPressure;
	delete[] middleDensity;
	delete[] middleVelocity;
	delete[] middlePressure;
}

void Simulation::initializeProfile(){
	downstreamR = 0;
	grid = new double[rgridNumber +1];
	pointDensityLeft = new double[rgridNumber + 1];
	pointDensityRight = new double[rgridNumber + 1];
	pointVelocity = new double[rgridNumber + 1];
	pointPressure = new double[rgridNumber + 1];
	middleDensity = new double[rgridNumber];
	middleVelocity = new double[rgridNumber];
	middlePressure = new double[rgridNumber];
	deltaR = (upstreamR - downstreamR)/rgridNumber;
	double r = downstreamR;
	double pressure0 = density0*kBoltzman*temperature/massProton;
	for(int i = 0; i < rgridNumber + 1; ++i){
		grid[i] = r;
		r += deltaR;
		switch(simulationType){
		case 1 :
			middleDensity[i] = density0;
			if(i > rgridNumber/100 && i < rgridNumber/10){
				middleVelocity[i] = U0;
			} else {
				middleVelocity[i] = 0;
			}
			middlePressure[i] = pressure0;
			break;
		case 2 :
			middleDensity[i] = density0 + density0*0.1*sin(i*10*2*pi/rgridNumber);
			middleVelocity[i] = sqrt(gamma*pressure0/density0)*density0*0.1*sin(i*10*2*pi/rgridNumber)/density0;
			middlePressure[i] = pressure0 + (gamma*pressure0/density0)*density0*0.1*sin(i*10*2*pi/rgridNumber);
			break;
		default:
			middleDensity[i] = density0;
			middleVelocity[i] = 0;
			if(i <= rgridNumber/100){
				middlePressure[i] = 100000000*pressure0;
			} else {
				middlePressure[i] = pressure0;
			}
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		pointDensityLeft[i] = middleDensity[i];
		pointDensityRight[i] = middleDensity[i];
		pointVelocity[i] = middleVelocity[i];
		pointPressure[i] = middlePressure[i];
	}

	if(cycleBound){
		pointDensityLeft[rgridNumber] = pointDensityLeft[0];
		pointDensityRight[rgridNumber] = pointDensityRight[0];
		pointVelocity[rgridNumber] = pointVelocity[0];
		pointPressure[rgridNumber] = pointPressure[0];
	} else {
		pointDensityLeft[rgridNumber] = pointDensityLeft[rgridNumber - 1];
		pointDensityRight[rgridNumber] = pointDensityRight[rgridNumber - 1];
		pointVelocity[rgridNumber] = pointVelocity[rgridNumber - 1];
		pointPressure[rgridNumber] = pointPressure[rgridNumber - 1];
	}
}

void Simulation::simulate(){
	FILE* outFile = fopen("./output/tamc_radial_profile.dat","w");
	FILE* outIteration = fopen("./output/iterations.dat","w");
	fclose(outIteration);
	printf("initialization\n");
	initializeProfile();
	output(outFile,this);
	fclose(outFile);
	updateMaxSoundSpeed();
	for(int i = 0; i < iterationNumber; ++i){
		printf("iteration ¹ %d\n", i);
		printf("time = %lf\n", time);
		printf("solving\n");
		time = time + deltaT;
		evaluateHydrodynamic();
		//updateValues();
		updateMaxSoundSpeed();
		//updateParameters();
		if(i % 1000 == 0){
			printf("outputing\n");
			outFile = fopen("./output/tamc_radial_profile.dat","a");
			output(outFile, this);
			fclose(outFile);
			outIteration = fopen("./output/iterations.dat","a");
			fprintf(outIteration, "%d %lf %lf\n", i, time, mass);
			fclose(outIteration);
		}
	}
}

void Simulation::evaluateHydrodynamic(){
	solveDiscontinious();

	for(int i = 0; i < rgridNumber; ++i){
		middleDensity[i] -= deltaT*(pointDensityLeft[i+1]*pointVelocity[i+1] - pointDensityRight[i]*pointVelocity[i])/deltaR;
		double middleMomentum = momentum(i) - deltaT*(pointPressure[i+1] + pointDensityLeft[i+1]*pointVelocity[i+1]*pointVelocity[i+1] - pointPressure[i] - pointDensityRight[i]*pointVelocity[i]*pointVelocity[i])/deltaR;
		double middleEnergy = energy(i) - deltaT*(pointPressure[i+1]*pointVelocity[i+1]*gamma/(gamma - 1) + pointDensityLeft[i+1]*pointVelocity[i+1]*pointVelocity[i+1]*pointVelocity[i+1]/2 - pointPressure[i]*pointVelocity[i]*gamma/(gamma - 1) - pointDensityRight[i]*pointVelocity[i]*pointVelocity[i]*pointVelocity[i]/2)/deltaR;
		middleVelocity[i] = middleMomentum/middleDensity[i];
		middlePressure[i] = (middleEnergy - middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2)*(gamma - 1);

		if(middleDensity[i] <= epsilon*density0){
			middleDensity[i] = epsilon*density0;
			middleVelocity[i] = 0;
			middlePressure[i] = 0;
		}
	}
}

void Simulation::solveDiscontinious(){

}

void Simulation::TracPen(double* u, double* flux, double cs, double leftFlux, double rightFlux){

}

double Simulation::momentum(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i < rgridNumber) {
		return middleDensity[i]*middleVelocity[i];
	} else {
		printf("i >= rgridNumber");
	}
	return 0;
}

double Simulation::energy(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i < rgridNumber) {
		return middlePressure[i]*(gamma - 1) + middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2;
	} else {
		printf("i >= rgridNumber");
	}
	return 0;
}

double Simulation::temperatureIn(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i < rgridNumber) {
		return middlePressure[i]*massProton/(kBoltzman*middleDensity[i]);
	} else {
		printf("i >= rgridNumber");
	}
	return 0;
}

double Simulation::minmod(double a, double b){
	return 0;
	if(a*b > 0){
		if(abs(a) < abs(b)){
			return a;
		} else {
			return b;
		}
	} else {
		return 0;
	}
}

void Simulation::updateMaxSoundSpeed(){
	maxSoundSpeed = 0;
	for(int i = 0; i < rgridNumber; ++i){
		double cs = (sqrt(gamma*middlePressure[i]/middleDensity[i]) + abs(middleVelocity[i]));
		if(cs > maxSoundSpeed){
			maxSoundSpeed = cs;
		} 
	}
	deltaT = 0.005*deltaR/maxSoundSpeed;
}
