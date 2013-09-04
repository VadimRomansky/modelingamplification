#include "stdafx.h"
#include "SpaceBin.h"
#include "particle.h"
#include "constants.h"
#include "simulation.h"
#include "util.h"
#include "matrix3d.h"


SpaceBin::SpaceBin(){
}

SpaceBin::SpaceBin(double R, double deltaR, double u, double rho, double T, double B, int i, bool scattering){
	r = R;
	r1 = r - deltaR/2;
	r2 = r + deltaR/2;
	U = u;
	density = rho;
	pressure = T*density*kBoltzman/massProton;
	temperature = T;
	B0 = B;
	numberR = i;
	smallAngleScattering = scattering;
	volume = 4*pi*(cube(r2) - cube(r1))/3;
}

SpaceBin::~SpaceBin(){
}

int SpaceBin::propagateParticle(Particle* particle ,double& time, double timeStep, const int rgridNumber){
	return 0;
}


int SpaceBin::binByCoordinates(double r, double theta, double phi, double r0, double deltar, double deltatheta, double deltaphi, int rgridNumber){
	int i = lowerInt((r - r0)/deltar);
	int j = lowerInt(theta/deltatheta);
	if( theta == pi ){
		theta = thetagridNumber - 1;
	}
	int k = lowerInt(phi/deltaphi);
	if(k >= phigridNumber){
		k = 0;
	}
	if(k < 0){
		k = phigridNumber - 1;
	}
	if (phi == 2*pi){
		k = 0;
	}
	if(i > rgridNumber - 1){
		//printf("aaa");
	}
	if( i < 0){
		//printf("aaa");
	}
	if(j >= thetagridNumber){
		j = 0;
	}
	if(j < 0){
		j = thetagridNumber - 1;
	}
	if(k > phigridNumber - 1){
		printf("aaa");
	}
	return i;
}

double SpaceBin::getFreePath(Particle* particle){
	//return speed_of_light*particle.localMomentum/(particle.Z*electron_charge*B);
	double lambda = speed_of_light*particle->localMomentum/(particle->Z*electron_charge*B0);
	if( lambda != lambda){
		printf("aaa");
	}
	if( 0*lambda != 0*lambda){
		printf("aaa");
	}
	return lambda;
}

void SpaceBin::makeOneStep(Particle* particle, double colisionTime, double& time){

}

void SpaceBin::multiplyParticleWeights(double value){

}

double SpaceBin::getEnergy(){
	return density*U*U/2 + pressure/(gamma - 1);
}


	


