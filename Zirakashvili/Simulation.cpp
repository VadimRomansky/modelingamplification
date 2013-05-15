#include "stdafx.h"
#include "simulation.h"
#include "util.h"
#include "constants.h"

Simulation::Simulation(){
	timeStep = defaultTimeStep;
	tau = defaultTimeStep;
}

void Simulation::initializeProfile(){
}

void Simulation::simulate(){
}

void Simulation::solveUpstream1(){
}

void Simulation::solveUpstream2(){
	double* density = new double[rgridNumber];
	double* velocity = new double[rgridNumber];
	double* pressure = new double[rgridNumber];
	for(int i = rgridNumber - 1; i >= 0; --i){
		double G = upstreamBins2[i]->density*exp(-3*tau) - density[i+1]*(3*tau*upstreamBins2[i+1]->xi*upstreamBins2[i+1]->xi*(velocity[i+1] - forwardV*upstreamBins2[i+1]->xi)/(qube(upstreamBins2[i+1]->xi) - qube(upstreamBins2[i]->xi)));
		velocity[i] = (upstreamBins2[i]->density*upstreamBins2[i]->U*exp(-3*tau) - density[i+1]*velocity[i+1]*3*tau*upstreamBins2[i+1]->xi*upstreamBins2[i+1]->xi/(forwardV*(qube(upstreamBins2[i+1]->xi) - qube(upstreamBins2[i]->xi)))-(pressure[i+1] - upstreamBins2[i+1]->pressure)*tau/(forwardV*(upstreamBins2[i+1]->xi - upstreamBins2[i]->xi)))/G;
		density[i] = G/(1 - 3*tau*upstreamBins2[i+1]->xi*upstreamBins2[i+1]->xi*(velocity[i] - forwardV*upstreamBins2[i]->xi)/(forwardV*(qube(upstreamBins2[i+1]->xi) - qube(upstreamBins2[i]->xi))));
		pressure[i] = (upstreamBins2[i]->pressure*exp(-3*tau) - pressure[i+1]*3*tau*sqr(upstreamBins2[i+1]->xi)*(velocity[i+1] - forwardV*upstreamBins2[i+1]->xi)/(forwardV*(qube(upstreamBins2[i+1]->xi) - qube(upstreamBins2[i]->xi))))
			(1 - 3*tau*(sqr(upstreamBins2[i]->xi)*(velocity[i] - forwardV*upstreamBins2[i]->xi) + (gamma - 1)*(sqr(upstreamBins2[i+1]->xi)*upstreamBins2[i+1]->U - sqr(upstreamBins2[i]->xi)*upstreamBins2[i]->U))/(forwardV*(qube(upstreamBins2[i+1]->xi) - qube(upstreamBins2[i]->xi))));
	}

	for(int i = 0; i < rgridNumber; ++i){
		upstreamBins2[i]->U = velocity[i];
		upstreamBins2[i]->density = density[i];
		upstreamBins2[i]->pressure = pressure[i];
	}

	delete[] density;
	delete[] velocity;
	delete[] pressure;
}

void Simulation::solveDownstream1(){
}

void Simulation::solveDownstream1(){
}
