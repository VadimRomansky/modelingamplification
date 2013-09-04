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
	deltaR = (downstreamR - upstreamR)/(rgridNumber );
	double R = upstreamR + deltaR/2;
	bins = new SpaceBin*[rgridNumber];
    for(int i = 0; i < rgridNumber; ++i){
		double density = density0;
		double u;
		if(i < shockWavePoint){
			u = U0;
			} else {
			u = 0;
		}
		bins[i] = new SpaceBin(R, deltaR, u, density, temperature,B0, i,smallAngleScattering);
		R = R + deltaR;
	}
}

void Simulation::simulate(){
	FILE* outFile = fopen("./output/zprofile.dat","w");
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
		printf("solving upstream\n");
		evaluateHydrodynamic();
		updateValues();
		updateMaxSoundSpeed();
		updateParameters();
		if(i % 1000 == 0){
			printf("outputing\n");
			outFile = fopen("./output/zprofile.dat","a");
			output(outFile, this);
			fclose(outFile);
			outIteration = fopen("./output/iterations.dat","a");
			fprintf(outIteration, "%d %lf %lf\n", i, time, mass);
			fclose(outIteration);
		}
	}
}

void Simulation::evaluateHydrodynamic(){	
	double* u1 = new double[rgridNumber];
	double* u2 = new double[rgridNumber];
	double* u3 = new double[rgridNumber];

	double* F1 = new double[rgridNumber];
	double* F2 = new double[rgridNumber];
	double* F3 = new double[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		u1[i] = bins[i]->density*bins[i]->volume;
		u2[i] = u1[i]*bins[i]->U;
		u3[i] = bins[i]->getEnergy()*bins[i]->volume;

		F1[i] = bins[i]->density*bins[i]->r*bins[i]->r*bins[i]->U;
		//todo pressure?
		F2[i] = F1[i]*bins[i]->U + bins[i]->pressure;
		F3[i] = bins[i]->U*bins[i]->getEnergy()*bins[i]->volume + bins[i]->U*bins[i]->pressure;
	}

	double leftFlux1 = 0;
	double leftFlux2 = 0;
	double leftFlux3 = 0;

	double rightFlux1 = 0;
	double rightFlux2 = 0;
	double rightFlux3 = 0;

	TracPen(u1, F1, maxSoundSpeed, leftFlux1, rightFlux1);
	TracPen(u2, F2, maxSoundSpeed, leftFlux2, rightFlux2);
	/*for(int i = 0; i < rgridNumber; ++i){
		u2[i] += tau*bins[i]->r*forwardShockWaveR*2*bins[i]->pressure*deltaF/forwardV;
	}*/
	TracPen(u3, F3, maxSoundSpeed, leftFlux3, rightFlux3);

	for(int i = 0; i < rgridNumber; ++i){
		bins[i]->u1 = u1[i];
		bins[i]->u2 = u2[i];
		bins[i]->u3 = u3[i];
	}

	delete[] u1;
	delete[] u2;
	delete[] u3;
	delete[] F1;
	delete[] F2;
	delete[] F3;

}

void Simulation::TracPen(double* u, double* flux, double cs, double leftFlux, double rightFlux){

	/*u[0] -= tau*((flux[0] + flux[1])/2 - leftFlux)/deltaXi;
	for(int i = 1; i < rgridNumber - 1; ++i){
		u[i] -= tau*0.5*(flux[i+1] - flux[i-1])/deltaXi;
	}
	u[rgridNumber - 1] -= tau*(rightFlux - (flux[rgridNumber - 1] + flux[rgridNumber - 2])/2)/deltaXi;*/
	tau = deltaT;
	double* uplus = new double[rgridNumber];
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

	u[0] -= tau*(0.5*(fplus[1] + fminus[1]) - leftFlux);
	for(int i = 1; i <= rgridNumber-2; ++i){
		u[i] -= tau*0.5*(fplus[i+1] + fminus[i+1] - (fplus[i] + fminus[i]));
	}
	u[rgridNumber - 1] -= tau*(rightFlux - 0.5*(fplus[rgridNumber - 1] + fminus[rgridNumber - 1]));

	delete[] uplus;
	delete[] uminus;

	delete[] fplus;
	delete[] fminus;
}

double Simulation::minmod(double a, double b){
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
	deltaT = 05*deltaR/maxSoundSpeed;
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
	}
}
