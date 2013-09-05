#include "stdafx.h"
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

Simulation::Simulation(){
	time = 0;
	timeStep = defaultTimeStep;
	tau = defaultTimeStep;
}

Simulation::~Simulation(){
	delete[] bins;
}

void Simulation::initializeProfile(){
	downstreamR = 0;
	deltaR = (upstreamR - downstreamR)/(rgridNumber - 1);
	double R = downstreamR + deltaR/2;
	bins = new SpaceBin*[rgridNumber];
    for(int i = 0; i < rgridNumber; ++i){
		double density = density0;
		switch(simulationType){
		case 1:
			//int shockWavePoint = rgridNumber/10;
			double u;
			if(i < rgridNumber/10){
				u = U0;
			} else {
				u = 0;
			}
			bins[i] = new SpaceBin(R, deltaR, u, density, temperature, B0, i, smallAngleScattering);
			R = R + deltaR;
			break;
		default:
			if(i == 0){
				bins[i] = new SpaceBin(R, deltaR, 0, density, 100000000, B0, i, smallAngleScattering);
			} else {
				bins[i] = new SpaceBin(R, deltaR, 0, density, temperature, B0, i, smallAngleScattering);
			}
			R = R + deltaR;
		}
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
		updateParameters();
		if(i % 30 == 0){
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
	double* newDensity = new double[rgridNumber];
	double* newVelocity = new double[rgridNumber];
	double* newPressure = new double[rgridNumber];
	
	newDensity[0] = bins[0]->density - deltaT*(2*bins[0]->U*bins[0]->density/bins[0]->r + bins[0]->density*(bins[1]->U - bins[0]->U)/deltaR + bins[0]->U*(bins[1]->density - bins[0]->density)/deltaR);
    newVelocity[0] = bins[0]->U - deltaT*(bins[0]->U*(bins[1]->U - bins[0]->U)/deltaR + (bins[1]->pressure - bins[0]->pressure)/(bins[0]->density*deltaR));
	newPressure[0] = bins[0]->pressure - deltaT*(bins[0]->U*(bins[1]->pressure - bins[0]->pressure)/deltaR + 2*gamma*bins[0]->pressure*bins[0]->U/bins[0]->r + gamma*bins[0]->pressure*(bins[1]->U - bins[0]->U)/deltaR);
	for(int i = 1; i < rgridNumber; ++i){
		newDensity[i] = bins[i]->density - deltaT*(2*bins[i]->U*bins[i]->density/bins[i]->r + bins[i]->density*(bins[i]->U - bins[i-1]->U)/deltaR + bins[i]->U*(bins[i]->density - bins[i-1]->density)/deltaR);
		newVelocity[i] = bins[i]->U - deltaT*(bins[i]->U*(bins[i]->U - bins[i-1]->U)/deltaR + (bins[i]->pressure - bins[i-1]->pressure)/(bins[i]->density*deltaR));
		newPressure[i] = bins[i]->pressure - deltaT*(bins[i]->U*(bins[i]->pressure - bins[i-1]->pressure)/deltaR + 2*gamma*bins[i]->pressure*bins[i]->U/bins[i]->r + gamma*bins[i]->pressure*(bins[i]->U - bins[i-1]->U)/deltaR);
	}

	for(int i = 0; i < rgridNumber; ++i){
		bins[i]->density = newDensity[i];
		bins[i]->U = newVelocity[i];
		bins[i]->pressure = newPressure[i];
	}

	delete[] newDensity;
	delete[] newPressure;
	delete[] newVelocity;
}

void Simulation::TracPen(double* u, double* flux, double cs, double leftFlux, double rightFlux){
	tau = deltaT;
	u[0] -= tau*((flux[0] + flux[1])/2 - leftFlux);
	for(int i = 1; i < rgridNumber - 1; ++i){
		u[i] -= tau*0.5*(flux[i] - flux[i-1]);
	}
	u[rgridNumber - 1] -= tau*(rightFlux - (flux[rgridNumber - 1] + flux[rgridNumber - 2])/2);
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
		double cs = (sqrt(gamma*bins[i]->pressure/bins[i]->density) + abs(bins[i]->U));
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
