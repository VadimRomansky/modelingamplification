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
	delete[] pointDensity;
	delete[] pointVelocity;
	delete[] pointPressure;
	delete[] middleDensity;
	delete[] middleMomentum;
	delete[] middleEnergy;
}

void Simulation::initializeProfile(){
	downstreamR = 0;
	grid = new double[rgridNumber +1];
	pointDensity = new double[rgridNumber + 1];
	pointVelocity = new double[rgridNumber + 1];
	pointPressure = new double[rgridNumber + 1];
	middleDensity = new double[rgridNumber];
	middleMomentum = new double[rgridNumber];
	middleEnergy = new double[rgridNumber];
	deltaR = (upstreamR - downstreamR)/rgridNumber;
	double r = downstreamR;
	double pressure0 = density0*kBoltzman*temperature/massProton;
	for(int i = 0; i < rgridNumber + 1; ++i){
		grid[i] = r;
		r += deltaR;
		switch(simulationType){
		case 1 :
			pointDensity[i] = density0;
			if(i > rgridNumber/100 && i < rgridNumber/10){
				pointVelocity[i] = U0;
			} else {
				pointVelocity[i] = 0;
			}
			pointPressure[i] = pressure0;
			break;
		case 2 :
			pointDensity[i] = density0 + density0*0.1*sin(i*10*2*pi/rgridNumber);
			pointVelocity[i] = sqrt(gamma*pressure0/density0)*density0*0.1*sin(i*10*2*pi/rgridNumber)/density0;
			pointPressure[i] = pressure0 + (gamma*pressure0/density0)*density0*0.1*sin(i*10*2*pi/rgridNumber);
			break;
		default:
			pointDensity[i] = density0;
			pointVelocity[i] = 0;
			if(i <= rgridNumber/100){
				pointPressure[i] = 100000000*pressure0;
			} else {
				pointPressure[i] = pressure0;
			}
		}
	}
	//cycle bounds
	if(cycleBound){
		pointDensity[rgridNumber] = pointDensity[0];
		pointVelocity[rgridNumber] = pointVelocity[0];
		pointPressure[rgridNumber] = pointPressure[0];
	}

	for(int i = 0; i < rgridNumber; ++i){
		middleDensity[i] = pointDensity[i+1];
		middleMomentum[i] = momentumAtPoint(i+1);
		middleEnergy[i] = energyAtPoint(i+1);
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
	for(int i = 0; i < rgridNumber; ++i){
		middleDensity[i] -= deltaT*(pointDensity[i+1]*pointVelocity[i+1] - pointDensity[i]*pointVelocity[i])/deltaR;
		middleMomentum[i] -= deltaT*(pointPressure[i+1] + pointDensity[i+1]*pointVelocity[i+1]*pointVelocity[i+1] - pointPressure[i] - pointDensity[i]*pointVelocity[i]*pointVelocity[i])/deltaR;
		middleEnergy[i] -= deltaT*(pointPressure[i+1]*pointVelocity[i+1]*gamma/(gamma - 1) + pointDensity[i+1]*pointVelocity[i+1]*pointVelocity[i+1]*pointVelocity[i+1]/2 - pointPressure[i]*pointVelocity[i]*gamma/(gamma - 1) - pointDensity[i]*pointVelocity[i]*pointVelocity[i]*pointVelocity[i]/2)/deltaR;
		if(middleDensity[i] <= epsilon*density0){
			middleDensity[i] = epsilon*density0;
			middleMomentum[i]= 0;
			middleEnergy[i] = 0;
		}
	}

	//at zero pint is constant
	for(int i = 1; i < rgridNumber+1; ++i){
		/*pointDensity[i] = 2*middleDensity[i-1] - pointDensity[i-1];
		double pointMomentum = 2*middleMomentum[i-1] - momentumAtPoint(i);
		double pointEnergy = 2*middleEnergy[i-1] - energyAtPoint(i);
		pointVelocity[i] = pointMomentum/pointDensity[i];
		pointPressure[i] = (pointEnergy - pointDensity[i]*pointVelocity[i]*pointVelocity[i]/2)*(gamma - 1)/gamma;*/
		pointDensity[i] = middleDensity[i-1];
		double pointMomentum = middleMomentum[i-1];
		double pointEnergy = middleEnergy[i-1];
		pointVelocity[i] = pointMomentum/pointDensity[i];
		pointPressure[i] = (pointEnergy - pointDensity[i]*pointVelocity[i]*pointVelocity[i]/2)*(gamma - 1)/gamma;
	}

	//cycle bond condition
	if(cycleBound){
		pointDensity[0] = pointDensity[rgridNumber];
		pointVelocity[0] = pointVelocity[rgridNumber];
		pointPressure[0] = pointPressure[rgridNumber];
	}
}

void Simulation::TracPen(double* u, double* flux, double cs, double leftFlux, double rightFlux){
	/*tau = deltaT;
	u[0] -= tau*((flux[0] + flux[1])/2 - leftFlux);
	for(int i = 1; i < rgridNumber - 1; ++i){
		u[i] -= tau*0.5*(flux[i] - flux[i-1]);
	}
	u[rgridNumber - 1] -= tau*(rightFlux - (flux[rgridNumber - 1] + flux[rgridNumber - 2])/2);*/
	/*double* uplus = new double[rgridNumber];
	double* uminus = new double[rgridNumber];

	double* fplus = new double[rgridNumber];
	double* fminus = new double[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		uplus[i] = cs*u[i] + flux[i];
		uminus[i] = cs*u[i] - flux[i];
	}

	fplus[0] = uplus[0];
	fminus[0] = -uminus[0];
	for(int i = 1; i < rgridNumber-1; ++i){

		if(i == 1){
			fplus[i] = uplus[i-1];
		} else {
			fplus[i] = uplus[i-1] + 0.5*minmod(uplus[i] - uplus[i-1], uplus[i-1] - uplus[i-2]);
		}

		if(i == rgridNumber - 1){
			fminus[i] = -uminus[i];
		} else {
		    fminus[i] = -uminus[i] + 0.5*minmod(uminus[i] - uminus[i-1], uminus[i+1] - uminus[i]);
		}

	}
	fplus[rgridNumber - 1] = uplus[rgridNumber - 2] + 0.5*minmod(uplus[rgridNumber - 1] - uplus[rgridNumber - 2], uplus[rgridNumber - 2] - uplus[rgridNumber - 3]);
	fminus[rgridNumber - 1] = -uminus[rgridNumber - 1];

	//u[0] -= tau*(0.5*(fplus[1] + fminus[1]) - leftFlux);
	u[0] -= tau*((flux[0] + flux[1])/2 - leftFlux);
	for(int i = 1; i <= rgridNumber-2; ++i){
		u[i] -= tau*0.5*(fplus[i+1] + fminus[i+1] - (fplus[i] + fminus[i]));
	}
	//u[rgridNumber - 1] -= tau*(rightFlux - 0.5*(fplus[rgridNumber - 1] + fminus[rgridNumber - 1]));
	u[rgridNumber - 1] -= tau*(rightFlux - (flux[rgridNumber - 1] + flux[rgridNumber - 2])/2);

	delete[] uplus;
	delete[] uminus;

	delete[] fplus;
	delete[] fminus;*/
}

double Simulation::momentumAtPoint(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i <= rgridNumber) {
		return pointDensity[i]*pointVelocity[i];
	} else if(i > rgridNumber) {
		printf("i > rgridNumber");
	}
	return 0;
}

double Simulation::energyAtPoint(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i <= rgridNumber) {
		return pointPressure[i]*gamma/(gamma - 1) + pointDensity[i]*pointVelocity[i]*pointVelocity[i]/2;
	} else if(i > rgridNumber) {
		printf("i > rgridNumber");
	}
	return 0;
}

double Simulation::temperatureAtPoint(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i <= rgridNumber) {
		return pointPressure[i]*massProton/(kBoltzman*pointDensity[i]);
	} else if(i > rgridNumber) {
		printf("i > rgridNumber");
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
		double cs = (sqrt(gamma*pointPressure[i]/pointDensity[i]) + abs(pointVelocity[i]));
		if(cs > maxSoundSpeed){
			maxSoundSpeed = cs;
		} 
	}
	deltaT = 0.005*deltaR/maxSoundSpeed;
}

void Simulation::updateParameters(){
	mass = 0;
	for(int i = 0; i < rgridNumber; ++i){
		bins[i]->temperature = bins[i]->pressure*massProton/(bins[i]->density*kBoltzman);
		mass += bins[i]->density*bins[i]->volume;
	}
}

void Simulation::updateValues(){
	for(int i = 0; i < rgridNumber; ++i){
		bins[i]->density = bins[i]->u1/(bins[i]->volume);
		bins[i]->U = bins[i]->u2/(bins[i]->density*bins[i]->volume);
		bins[i]->pressure = (gamma - 1)*(bins[i]->u3/(bins[i]->volume) - bins[i]->density*bins[i]->U*bins[i]->U/2);
		if(bins[i]->density < 0){
			printf("density < 0 bin number %d\n", i);
			bins[i]->density = density0*epsilon;
			bins[i]->U = 0;
			bins[i]->pressure = 0;
		}
	}
}
