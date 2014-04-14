#include "simulation.h"
#include "constants.h"
#include "util.h"

void Simulation::evauateField(){

}

void Simulation::evaluateCRFlux(double* crflux){
}

void Simulation::growthRate(double* crflux, double* rate){
	double J = 0;
	for(int i = 0; i < pgridNumber; ++i){
		J = crflux[i]*cube(pgrid[i])*deltaLogP;
	}

	for(int i = 0; i < kgridNumber; ++i){
		double z = kgrid[i]*speed_of_light*pgrid[j]/(electron_charge*Bls);
	}
}