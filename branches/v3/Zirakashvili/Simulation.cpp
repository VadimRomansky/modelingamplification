#include "stdafx.h"
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

Simulation::Simulation(){
	timeStep = defaultTimeStep;
	tau = defaultTimeStep;
}

void Simulation::initializeProfile(){
	forwardV = U0;
	reverseV = 0.5*forwardV;
	contactDiscontV = 0.75*forwardV;

	forwardShockWaveR = 0.5*upstreamR;
	contactDiscontR = 0.95*forwardShockWaveR;
	reverseShockWaveR = 0.9*forwardShockWaveR;

	double p0 = density0*kBoltzman*temperature/massProton;
	for(int i = 0; i < rgridNumber; ++i){
		upstreamBins2[i]->r = forwardShockWaveR + (i + 0.5)*(upstreamR - forwardShockWaveR)/rgridNumber;
		upstreamBins2[i]->xi = upstreamBins2[i]->r/forwardShockWaveR;
		upstreamBins2[i]->density = density0;
		upstreamBins2[i]->U = 0;
		upstreamBins2[i]->pressure = p0;

		downstreamBins2[i]->r = contactDiscontR + (i + 0.5)*(forwardShockWaveR - contactDiscontR)/rgridNumber;
		downstreamBins2[i]->xi = (downstreamBins1[i]->r - contactDiscontR)/(forwardShockWaveR - contactDiscontR);
		downstreamBins2[i]->density = 4*density0;
		downstreamBins2[i]->U = 0.75*U0;
		downstreamBins2[i]->pressure = 0.75*density0*U0*U0;

		downstreamBins1[i]->r = reverseShockWaveR +(i + 0.5)*(contactDiscontR - reverseShockWaveR)/rgridNumber;
		downstreamBins1[i]->xi = (downstreamBins1[i]->r - contactDiscontR)/(contactDiscontR - reverseShockWaveR);
		downstreamBins1[i]->density = 4*density0;
		downstreamBins1[i]->U = 0.75*U0;
		downstreamBins1[i]->pressure = 0.75*density0*U0*U0;

		upstreamBins1[i]->r = reverseShockWaveR*(i + 0.5)/rgridNumber;
		upstreamBins1[i]->xi = upstreamBins1[i]->r/reverseShockWaveR;
		upstreamBins1[i]->density = density0*sqr(upstreamBins2[i]->r/reverseShockWaveR);
		upstreamBins1[i]->U = 1.5*forwardV*upstreamBins1[i]->r/reverseShockWaveR;
		upstreamBins1[i]->pressure = 0;
	}
}

void Simulation::simulate(){
	FILE* outFile = fopen("./output/zprofile.dat","w");
	printf("initialization\n");
	initializeProfile();
	output(outFile,this);
	fclose(outFile);
	moveShockWaves();
	for(int i = 0; i < iterationNumber; ++i){
		printf("iteration ¹ %d", i);
		printf("solving upstream\n");
		solveUpstream1();
		solveUpstream2();
		printf("solving downstream\n");
		solveDownstream1();
		solveDownstream2();
		printf("moving shock waves\n");
		moveShockWaves();
		printf("outputing\n");
		outFile = fopen("./output/zprofile.dat","a");
		output(outFile, this);
		fclose(outFile);
	}
}

void Simulation::solveUpstream1(){
	double* density = new double[rgridNumber];
	double* velocity = new double[rgridNumber];
	double* pressure = new double[rgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		double xi1 = upstreamBins1[i]->xi;
		double xi2 = upstreamBins1[i-1]->xi;
		double H1 = 1 + forwardShockWaveR*tau*(upstreamBins1[1]->U - reverseV*xi1)/(reverseShockWaveR*forwardV*(xi1 - xi2));
		double H2 = upstreamBins1[i]->density*cube(oldReverseShockWaveR/reverseShockWaveR) + density[i-1]*3*forwardShockWaveR*tau*xi2*xi2*(velocity[i-1]-reverseV*xi2)/(reverseShockWaveR*forwardV*(cube(xi1) - cube(xi2)));
		velocity[i] = (upstreamBins1[i]->U + velocity[i-1]*forwardShockWaveR*tau*(upstreamBins1[i]->U - reverseV*xi1)/(reverseShockWaveR*forwardV*(xi1 - xi2)) - forwardShockWaveR*tau*(upstreamBins1[i]->pressure - upstreamBins1[i-1]->pressure)/(reverseShockWaveR*forwardV*(xi1 - xi2)))/H1;
		//todo u old or new?
		density[i] = H2/(1 + 3*forwardShockWaveR*tau*xi1*xi1*(velocity[i] - reverseV*xi1)/(reverseShockWaveR*forwardV*(cube(xi1) - cube(xi2))));
		pressure[i] = (upstreamBins1[i]->pressure*cube(oldReverseShockWaveR/reverseShockWaveR) + pressure[i-1]*3*forwardShockWaveR*tau*xi2*xi2*(velocity[i-1] - reverseV*xi2)/(reverseShockWaveR*forwardV*(cube(xi1) - cube(xi2))))/
			(1 + 3*forwardShockWaveR*tau*(xi1*xi1*(velocity[i] - reverseV*xi1) + (gamma - 1)*(xi1*xi1*velocity[i] - xi2*xi2*velocity[i-1]))/(reverseShockWaveR*forwardV*(cube(xi1) - cube(xi2))));
	}

	for(int i = 0; i < rgridNumber; ++i){
		upstreamBins1[i]->U = velocity[i];
		upstreamBins1[i]->density = density[i];
		upstreamBins1[i]->pressure = pressure[i];
	}

	delete[] density;
	delete[] velocity;
	delete[] pressure;
}

void Simulation::solveUpstream2(){
	double* density = new double[rgridNumber];
	double* velocity = new double[rgridNumber];
	double* pressure = new double[rgridNumber];
	for(int i = rgridNumber - 1; i >= 0; --i){
		double G = upstreamBins2[i]->density*cube(oldForwardShockWaveR/forwardShockWaveR) - density[i+1]*(3*tau*upstreamBins2[i+1]->xi*upstreamBins2[i+1]->xi*(velocity[i+1] - forwardV*upstreamBins2[i+1]->xi)/(forwardV*(cube(upstreamBins2[i+1]->xi) - cube(upstreamBins2[i]->xi))));
		velocity[i] = (upstreamBins2[i]->density*upstreamBins2[i]->U*cube(oldForwardShockWaveR/forwardShockWaveR) - density[i+1]*velocity[i+1]*3*tau*upstreamBins2[i+1]->xi*upstreamBins2[i+1]->xi*(velocity[i+1]-forwardV*upstreamBins2[i+1]->xi)/(forwardV*(cube(upstreamBins2[i+1]->xi) - cube(upstreamBins2[i]->xi)))-(upstreamBins2[i+1]->pressure - upstreamBins2[i]->pressure)*tau/(forwardV*(upstreamBins2[i+1]->xi - upstreamBins2[i]->xi)))/G;
		density[i] = G/(1 - 3*tau*upstreamBins2[i+1]->xi*upstreamBins2[i+1]->xi*(velocity[i] - forwardV*upstreamBins2[i]->xi)/(forwardV*(cube(upstreamBins2[i+1]->xi) - cube(upstreamBins2[i]->xi))));
		pressure[i] = (upstreamBins2[i]->pressure*cube(oldForwardShockWaveR/forwardShockWaveR) - pressure[i+1]*3*tau*sqr(upstreamBins2[i+1]->xi)*(velocity[i+1] - forwardV*upstreamBins2[i+1]->xi)/(forwardV*(cube(upstreamBins2[i+1]->xi) - cube(upstreamBins2[i]->xi))))*
			(1 - 3*tau*(sqr(upstreamBins2[i]->xi)*(velocity[i] - forwardV*upstreamBins2[i]->xi) + (gamma - 1)*(sqr(upstreamBins2[i+1]->xi)*velocity[i+1] - sqr(upstreamBins2[i]->xi)*velocity[i]))/(forwardV*(cube(upstreamBins2[i+1]->xi) - cube(upstreamBins2[i]->xi))));
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
	double deltaB = contactDiscontR - reverseShockWaveR;
	double deltaXi = deltaB/rgridNumber;
	double* u1 = new double[rgridNumber];
	double* u2 = new double[rgridNumber];
	double* u3 = new double[rgridNumber];

	double* F1 = new double[rgridNumber];
	double* F2 = new double[rgridNumber];
	double* F3 = new double[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		double vb = downstreamBins1[i]->U - contactDiscontV*(1+downstreamBins1[i]->xi) + downstreamBins1[i]->xi*reverseV;
		F1[i] = forwardShockWaveR*downstreamBins1[i]->density*downstreamBins1[i]->r*downstreamBins1[i]->r*vb/forwardV;
		F2[i] = forwardShockWaveR*downstreamBins1[i]->r*downstreamBins1[i]->r*(downstreamBins1[i]->density*downstreamBins1[i]->U*vb + downstreamBins1[i]->pressure)/forwardV;
		F3[i] = forwardShockWaveR*downstreamBins1[i]->r*downstreamBins1[i]->r*(downstreamBins1[i]->getEnergy()*vb + downstreamBins1[i]->U*downstreamBins1[i]->pressure)/forwardV;
	}


	TracPen(u1, F1, maxSoundSpeed, deltaXi);
	TracPen(u2, F2, maxSoundSpeed, deltaXi);
	for(int i = 0; i < rgridNumber; ++i){
		u2[i] += tau*downstreamBins1[i]->r*forwardShockWaveR*2*downstreamBins1[i]->pressure*deltaB/forwardV;
	}
	TracPen(u3, F3, maxSoundSpeed, deltaXi);

	delete[] u1;
	delete[] u2;
	delete[] u3;
	delete[] F1;
	delete[] F2;
	delete[] F3;
}

void Simulation::solveDownstream2(){	double deltaB = contactDiscontR - reverseShockWaveR;
	double deltaF = forwardShockWaveR - contactDiscontR;
	double deltaXi = deltaF/rgridNumber;
	double* u1 = new double[rgridNumber];
	double* u2 = new double[rgridNumber];
	double* u3 = new double[rgridNumber];

	double* F1 = new double[rgridNumber];
	double* F2 = new double[rgridNumber];
	double* F3 = new double[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		double vf = downstreamBins2[i]->U - contactDiscontV*(1-downstreamBins2[i]->xi) - downstreamBins2[i]->xi*forwardV;
		F1[i] = forwardShockWaveR*downstreamBins2[i]->density*downstreamBins2[i]->r*downstreamBins2[i]->r*vf/forwardV;
		F2[i] = forwardShockWaveR*downstreamBins2[i]->r*downstreamBins2[i]->r*(downstreamBins2[i]->density*downstreamBins2[i]->U*vf + downstreamBins2[i]->pressure)/forwardV;
		F3[i] = forwardShockWaveR*downstreamBins2[i]->r*downstreamBins2[i]->r*(downstreamBins2[i]->getEnergy()*vf + downstreamBins2[i]->U*downstreamBins2[i]->pressure)/forwardV;
	}

	TracPen(u1, F1, maxSoundSpeed, deltaXi);
	TracPen(u2, F2, maxSoundSpeed, deltaXi);
	for(int i = 0; i < rgridNumber; ++i){
		u2[i] += tau*downstreamBins1[i]->r*forwardShockWaveR*2*downstreamBins2[i]->pressure*deltaF/forwardV;
	}
	TracPen(u3, F3, maxSoundSpeed, deltaXi);

	delete[] u1;
	delete[] u2;
	delete[] u3;
	delete[] F1;
	delete[] F2;
	delete[] F3;

}

void Simulation::TracPen(double* u, double* flux, double cs, double deltaXi){
	double* uplus = new double[rgridNumber];
	double* uminus = new double[rgridNumber];

	double* fplus = new double[rgridNumber];
	double* fminus = new double[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		uplus[i] = cs*u[i] + flux[i];
		uminus[i] = cs*u[i] - flux[i];

		fplus[i] = uplus[i-1] + 0.5*minmod(uplus[i] - uplus[i-1], uplus[i-1] - uplus[i-2]);
		fminus[i] = -uplus[i] + 0.5*minmod(uminus[i] - uminus[i-1], uminus[i+1] - uminus[i]);

	}

	for(int i = 0; i < rgridNumber; ++i){
		u[i] -= tau*0.5*(fplus[i+1] + fminus[i+1] - fplus[i] - fminus[i])/deltaXi;;
	}

	delete[] uplus;
	delete[] uminus;

	delete[] fplus;
	delete[] uminus;
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

void Simulation::moveShockWaves(){
	double oldForwardV = forwardV;
	double oldReverseV = reverseV;
	double oldContactDiscontV = contactDiscontV;

	oldForwardShockWaveR = forwardShockWaveR;
	oldReverseShockWaveR = reverseShockWaveR;
	oldContactDiscontR = contactDiscontR;

	forwardV = upstreamBins2[0]->U + sqrt((gamma + 1)*downstreamBins2[rgridNumber - 1]->pressure/2 + (gamma - 1)*upstreamBins2[0]->pressure/2)/sqrt(upstreamBins2[0]->density);
	reverseV = upstreamBins1[rgridNumber - 1]->U - sqrt((gamma + 1)*downstreamBins1[0]->pressure/2 + (gamma - 1)*upstreamBins1[rgridNumber - 1]->pressure)/sqrt(upstreamBins1[rgridNumber-1]->density);
	SpaceBin* leftBin = downstreamBins1[rgridNumber - 1];
	SpaceBin* rightBin = downstreamBins2[0];
	contactDiscontV = (leftBin->U*sqrt(leftBin->density) + rightBin->U*sqrt(rightBin->density) + (leftBin->pressure - rightBin->pressure)/(sqrt(gamma*(leftBin->pressure + rightBin->pressure)/2)))/(sqrt(leftBin->density) + sqrt(rightBin->density));

	forwardShockWaveR *= exp(tau);
	reverseShockWaveR += oldForwardShockWaveR*tau*(reverseV + oldReverseV)*exp(0.5*tau)/(forwardV + oldForwardV);
	time += 2*oldForwardShockWaveR*tau*exp(0.5*tau)/(forwardV + oldForwardV);
}

void Simulation::updateMaxSoundSpeed(){
	double deltaF = forwardShockWaveR - contactDiscontR;
	double deltaB = contactDiscontR - reverseShockWaveR;
	maxSoundSpeed = 0;
	for(int i = 0; i < rgridNumber; ++i){
		double cs = forwardShockWaveR*(sqrt(gamma*downstreamBins1[i]->pressure/downstreamBins1[i]->density) + abs(downstreamBins1[i]->U))/(forwardV*deltaF);
		if(cs > maxSoundSpeed){
			maxSoundSpeed = cs;
		} 
		cs = forwardShockWaveR*(sqrt(gamma*downstreamBins2[i]->pressure/downstreamBins2[i]->density) + abs(downstreamBins2[i]->U))/(forwardV*deltaB);
		if(cs > maxSoundSpeed){
			maxSoundSpeed = cs;
		} 
	}

	tau = 0.5*max(deltaF,deltaB)/(rgridNumber*maxSoundSpeed);
}
