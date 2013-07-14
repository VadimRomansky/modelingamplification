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
	delete[] upstreamBins2;
	delete[] downstreamBins2;
	delete[] downstreamBins1;
	delete[] upstreamBins1;
}

void Simulation::initializeProfile(){
	upstreamBins2 = new SpaceBin*[rgridNumber];
	downstreamBins2 = new SpaceBin*[rgridNumber];
	downstreamBins1 = new SpaceBin*[rgridNumber];
	upstreamBins1 = new SpaceBin*[rgridNumber];

	forwardV = U0;
	reverseV = 0.5*forwardV;
	contactDiscontV = 0.75*forwardV;

	forwardShockWaveR = 0.5*upstreamR;
	contactDiscontR = 0.95*forwardShockWaveR;
	reverseShockWaveR = 0.9*forwardShockWaveR;

	double p0 = 2*(1 - (gamma + 1)*0.75/2)*density0*U0*U0/(gamma - 1);
	double xiMin = 0.000001;
	double xiMax = upstreamR/forwardShockWaveR;
	double k2 = log((xiMax + xiMin - 1)/xiMin)/(rgridNumber - 1);
	double k1 = 1.0/(rgridNumber - 1);

	int ejectaR = rgridNumber - 500;
	for(int i = 0; i < rgridNumber; ++i){

		upstreamBins2[i] = new SpaceBin();
		//upstreamBins2[i]->r = forwardShockWaveR + (i + 0.5)*(upstreamR - forwardShockWaveR)/rgridNumber;
		upstreamBins2[i]->xi = 1 - xiMin + xiMin*exp(i*k2);
		upstreamBins2[i]->r = upstreamBins2[i]->xi*forwardShockWaveR;
		upstreamBins2[i]->density = density0;
		upstreamBins2[i]->U = 0;
		upstreamBins2[i]->pressure = p0;
		//upstreamBins2[i]->temperature = upstreamBins2[i]->pressure*massProton/(upstreamBins2[i]->density*kBoltzman);
		upstreamBins2[i]->temperature = 10000;
		upstreamBins2[i]->pressure = upstreamBins2[i]->density*kBoltzman*upstreamBins2[i]->temperature/massProton;

		downstreamBins2[i] = new SpaceBin();
		downstreamBins2[i]->r = contactDiscontR + (i + 0.5)*(forwardShockWaveR - contactDiscontR)/rgridNumber;
		downstreamBins2[i]->xi = (downstreamBins2[i]->r - contactDiscontR)/(forwardShockWaveR - contactDiscontR);
		downstreamBins2[i]->density = 4*density0;
		downstreamBins2[i]->U = 0.75*U0;
		downstreamBins2[i]->pressure = 0.75*density0*U0*U0;
		downstreamBins2[i]->temperature = downstreamBins2[i]->pressure*massProton/(downstreamBins2[i]->density*kBoltzman);

		downstreamBins1[i] = new SpaceBin();
		downstreamBins1[i]->r = reverseShockWaveR + (i + 0.5)*(contactDiscontR - reverseShockWaveR)/rgridNumber;
		downstreamBins1[i]->xi = (downstreamBins1[i]->r - contactDiscontR)/(contactDiscontR - reverseShockWaveR);
		downstreamBins1[i]->density = 4*density0;
		downstreamBins1[i]->U = 0.75*U0;
		downstreamBins1[i]->pressure = 0.75*density0*U0*U0;
		downstreamBins1[i]->temperature = downstreamBins1[i]->pressure*massProton/(downstreamBins1[i]->density*kBoltzman);

		upstreamBins1[i] = new SpaceBin();
		//upstreamBins1[i]->r = reverseShockWaveR*(i + 0.5)/rgridNumber;
		upstreamBins1[i]->xi = 2*i*k1 + (log(1 + power(2*xiMin, 1 - 2*i*k1)) - log(1 + 2*xiMin))/log(2*xiMin);
		upstreamBins1[i]->r = upstreamBins1[i]->xi*reverseShockWaveR;
		//if( i > ejectaR){
		//	upstreamBins1[i]->density = density0*power(upstreamBins1[i]->r/reverseShockWaveR, -7);
		//} else {
		//}
		upstreamBins1[i]->density = density0;
		upstreamBins1[i]->U = 1.5*forwardV*upstreamBins1[i]->r/reverseShockWaveR;
		upstreamBins1[i]->pressure = 0;
		upstreamBins1[i]->temperature = upstreamBins1[i]->pressure*massProton/(upstreamBins1[i]->density*kBoltzman);
	}
	/*for(int i = 0; i <= ejectaR; ++i){
		upstreamBins1[i]->density = density0*power(upstreamBins1[ejectaR]->r/reverseShockWaveR, -7);
	}*/
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
	moveShockWaves();
	updateMaxSoundSpeed();
	for(int i = 0; i < iterationNumber; ++i){
		printf("iteration ¹ %d\n", i);
		printf("time = %lf\n", time);
		printf("solving upstream\n");
		solveUpstream1();
		solveUpstream2();
		printf("solving downstream\n");
		solveDownstream1();
		solveDownstream2();
		printf("moving shock waves\n");
		moveShockWaves();
		updateMaxSoundSpeed();
		updateDownstreamValues();
		updateParameters();
		if(i % 1000 == 0){
			printf("outputing\n");
			outFile = fopen("./output/zprofile.dat","a");
			output(outFile, this);
			fclose(outFile);
			outIteration = fopen("./output/iterations.dat","a");
			fprintf(outIteration, "%d %lf %lf %lf %lf %lf %lf\n", i, time, mass, upstreamMass1, downstreamMass1, downstreamMass2, upstreamMass2);
			fclose(outIteration);
		}
	}
}

void Simulation::solveUpstream1(){
	double* density = new double[rgridNumber];
	double* velocity = new double[rgridNumber];
	double* pressure = new double[rgridNumber];

	double dt = oldForwardShockWaveR*tau*exp(0.5*tau)/forwardV;

	//density[0] = upstreamBins1[0]->density;
	density[0] = upstreamBins1[0]->density - 3*dt*upstreamBins1[0]->density*upstreamBins1[1]->U/(upstreamBins1[1]->r);
	velocity[0] = upstreamBins1[0]->U;
	//pressure[0] = upstreamBins1[0]->pressure;
	pressure[0] = upstreamBins1[0]->pressure*density[0]/upstreamBins1[0]->density;
	for(int i = 1; i < rgridNumber; ++i){
		double xi1 = upstreamBins1[i]->xi;
		double xi2 = upstreamBins1[i-1]->xi;
		double H1 = 1 + forwardShockWaveR*tau*(upstreamBins1[i]->U - reverseV*xi1)/(reverseShockWaveR*forwardV*(xi1 - xi2));
		velocity[i] = (upstreamBins1[i]->U + velocity[i-1]*forwardShockWaveR*tau*(upstreamBins1[i]->U - reverseV*xi1)/(reverseShockWaveR*forwardV*(xi1 - xi2)) - forwardShockWaveR*tau*(upstreamBins1[i]->pressure - upstreamBins1[i-1]->pressure)/(reverseShockWaveR*forwardV*(xi1 - xi2)))/H1;

		double H2 = upstreamBins1[i]->density*cube(oldReverseShockWaveR/reverseShockWaveR) + density[i-1]*3*forwardShockWaveR*tau*xi2*xi2*(velocity[i-1]-reverseV*xi2)/(reverseShockWaveR*forwardV*(cube(xi1) - cube(xi2)));
		//todo u old or new?
		density[i] = H2/(1 + 3*forwardShockWaveR*tau*xi1*xi1*(upstreamBins1[i]->U - reverseV*xi1)/(reverseShockWaveR*forwardV*(cube(xi1) - cube(xi2))));
		if(density[i] <= 0){
			printf("density < 0 i = %d \n", i);
			density[i] = epsilon*density0;
		}
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
	density[rgridNumber - 1] = upstreamBins2[rgridNumber - 1]->density;
	velocity[rgridNumber - 1] = upstreamBins2[rgridNumber - 1]->U;
	pressure[rgridNumber - 1] = upstreamBins2[rgridNumber - 1]->pressure;
	for(int i = rgridNumber - 2; i >= 0; --i){
		double G = upstreamBins2[i]->density*cube(oldForwardShockWaveR/forwardShockWaveR) - density[i+1]*(3*tau*upstreamBins2[i+1]->xi*upstreamBins2[i+1]->xi*(velocity[i+1] - forwardV*upstreamBins2[i+1]->xi)/(forwardV*(cube(upstreamBins2[i+1]->xi) - cube(upstreamBins2[i]->xi))));
		velocity[i] = (upstreamBins2[i]->density*upstreamBins2[i]->U*cube(oldForwardShockWaveR/forwardShockWaveR) - density[i+1]*velocity[i+1]*3*tau*upstreamBins2[i+1]->xi*upstreamBins2[i+1]->xi*(velocity[i+1]-forwardV*upstreamBins2[i+1]->xi)/(forwardV*(cube(upstreamBins2[i+1]->xi) - cube(upstreamBins2[i]->xi)))-(upstreamBins2[i+1]->pressure - upstreamBins2[i]->pressure)*tau/(forwardV*(upstreamBins2[i+1]->xi - upstreamBins2[i]->xi)))/G;
		density[i] = G/(1 - 3*tau*upstreamBins2[i+1]->xi*upstreamBins2[i+1]->xi*(velocity[i] - forwardV*upstreamBins2[i]->xi)/(forwardV*(cube(upstreamBins2[i+1]->xi) - cube(upstreamBins2[i]->xi))));
		if(density[i] <= 0){
			printf("density < 0 i = %d \n", i);
			density[i] = epsilon*density0;
		}
		pressure[i] = (upstreamBins2[i]->pressure*cube(oldForwardShockWaveR/forwardShockWaveR) - pressure[i+1]*3*tau*sqr(upstreamBins2[i+1]->xi)*(velocity[i+1] - forwardV*upstreamBins2[i+1]->xi)/(forwardV*(cube(upstreamBins2[i+1]->xi) - cube(upstreamBins2[i]->xi))))/
			(1 - 3*tau*(sqr(upstreamBins2[i]->xi)*(velocity[i] - forwardV*upstreamBins2[i]->xi) + (gamma - 1)*(sqr(upstreamBins2[i+1]->xi)*velocity[i+1] - sqr(upstreamBins2[i]->xi)*velocity[i]))/(forwardV*(cube(upstreamBins2[i+1]->xi) - cube(upstreamBins2[i]->xi))));
		alertNaNOrInfinity(pressure[i],"pressure\n");
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
	double deltaXi = 1.0/(rgridNumber-1);
	double* u1 = new double[rgridNumber];
	double* u2 = new double[rgridNumber];
	double* u3 = new double[rgridNumber];

	double* F1 = new double[rgridNumber];
	double* F2 = new double[rgridNumber];
	double* F3 = new double[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		double vb = downstreamBins1[i]->U - contactDiscontV*(1+downstreamBins1[i]->xi) + downstreamBins1[i]->xi*reverseV;
		//double vb = downstreamBins1[i]->U;
		u1[i] = downstreamBins1[i]->density*downstreamBins1[i]->r*downstreamBins1[i]->r*deltaB;
		u2[i] = u1[i]*downstreamBins1[i]->U;
		u3[i] = downstreamBins1[i]->r*downstreamBins1[i]->r*downstreamBins1[i]->getEnergy()*deltaB;

		F1[i] = forwardShockWaveR*downstreamBins1[i]->density*downstreamBins1[i]->r*downstreamBins1[i]->r*vb/forwardV;
		F2[i] = forwardShockWaveR*downstreamBins1[i]->r*downstreamBins1[i]->r*(downstreamBins1[i]->density*downstreamBins1[i]->U*vb + downstreamBins1[i]->pressure)/forwardV;
		F3[i] = forwardShockWaveR*downstreamBins1[i]->r*downstreamBins1[i]->r*(downstreamBins1[i]->getEnergy()*vb + downstreamBins1[i]->U*downstreamBins1[i]->pressure)/forwardV;
	}

    double middlePressure = (downstreamBins1[rgridNumber - 1]->pressure*sqrt(downstreamBins2[0]->density) + downstreamBins2[0]->pressure*sqrt(downstreamBins1[rgridNumber - 1]->density) + (downstreamBins1[rgridNumber - 1]->U - downstreamBins2[0]->U)*sqrt(gamma*downstreamBins1[rgridNumber - 1]->density*downstreamBins2[0]->density*(downstreamBins1[rgridNumber - 1]->pressure + downstreamBins2[0]->pressure)/2))/(sqrt(downstreamBins1[rgridNumber - 1]->density) + sqrt(downstreamBins2[0]->density));

	double leftFlux1 = upstreamBins1[rgridNumber - 1]->density*sqr(reverseShockWaveR)*forwardShockWaveR*(upstreamBins1[rgridNumber - 1]->U - reverseV)/forwardV;
	double leftFlux2 = sqr(reverseShockWaveR)*forwardShockWaveR*(upstreamBins1[rgridNumber - 1]->density*upstreamBins1[rgridNumber - 1]->U*(upstreamBins1[rgridNumber - 1]->U - reverseV) + upstreamBins1[rgridNumber - 1]->pressure)/forwardV;
	double leftFlux3 = sqr(forwardShockWaveR)*forwardShockWaveR*(upstreamBins1[rgridNumber - 1]->getEnergy()*(upstreamBins1[rgridNumber - 1]->U - reverseV) + upstreamBins1[rgridNumber - 1]->U*upstreamBins1[rgridNumber - 1]->pressure)/forwardV;

	double rightFlux1 = 0;
	double rightFlux2 = middlePressure*contactDiscontR*contactDiscontR*forwardShockWaveR/forwardV;
	double rightFlux3 = contactDiscontV*rightFlux2;


	TracPen(u1, F1, maxSoundSpeed, deltaXi, leftFlux1, rightFlux1);
	TracPen(u2, F2, maxSoundSpeed, deltaXi, leftFlux2, rightFlux2);
	for(int i = 0; i < rgridNumber; ++i){
		u2[i] += tau*downstreamBins1[i]->r*forwardShockWaveR*2*downstreamBins1[i]->pressure*deltaB/forwardV;
	}
	TracPen(u3, F3, maxSoundSpeed, deltaXi, leftFlux3, rightFlux3);

	for(int i = 0; i < rgridNumber; ++i){
		downstreamBins1[i]->u1 = u1[i];
		downstreamBins1[i]->u2 = u2[i];
		downstreamBins1[i]->u3 = u3[i];
	}

	delete[] u1;
	delete[] u2;
	delete[] u3;
	delete[] F1;
	delete[] F2;
	delete[] F3;
}

void Simulation::solveDownstream2(){	
	double deltaF = forwardShockWaveR - contactDiscontR;
	double deltaXi = 1.0/(rgridNumber-1);
	double* u1 = new double[rgridNumber];
	double* u2 = new double[rgridNumber];
	double* u3 = new double[rgridNumber];

	double* F1 = new double[rgridNumber];
	double* F2 = new double[rgridNumber];
	double* F3 = new double[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		double vf = downstreamBins2[i]->U - contactDiscontV*(1-downstreamBins2[i]->xi) - downstreamBins2[i]->xi*forwardV;
		//double vf = downstreamBins2[i]->U;
		u1[i] = downstreamBins2[i]->density*downstreamBins2[i]->r*downstreamBins2[i]->r*deltaF;
		u2[i] = u1[i]*downstreamBins2[i]->U;
		u3[i] = downstreamBins2[i]->r*downstreamBins2[i]->r*downstreamBins2[i]->getEnergy()*deltaF;

		F1[i] = forwardShockWaveR*downstreamBins2[i]->density*downstreamBins2[i]->r*downstreamBins2[i]->r*vf/forwardV;
		F2[i] = forwardShockWaveR*downstreamBins2[i]->r*downstreamBins2[i]->r*(downstreamBins2[i]->density*downstreamBins2[i]->U*vf + downstreamBins2[i]->pressure)/forwardV;
		F3[i] = forwardShockWaveR*downstreamBins2[i]->r*downstreamBins2[i]->r*(downstreamBins2[i]->getEnergy()*vf + downstreamBins2[i]->U*downstreamBins2[i]->pressure)/forwardV;
	}

	double middlePressure = (downstreamBins1[rgridNumber - 1]->pressure*sqrt(downstreamBins2[0]->density) + downstreamBins2[0]->pressure*sqrt(downstreamBins1[rgridNumber - 1]->density) + (downstreamBins1[rgridNumber - 1]->U - downstreamBins2[0]->U)*sqrt(gamma*downstreamBins1[rgridNumber - 1]->density*downstreamBins2[0]->density*(downstreamBins1[rgridNumber - 1]->pressure + downstreamBins2[0]->pressure)/2))/(sqrt(downstreamBins1[rgridNumber - 1]->density) + sqrt(downstreamBins2[0]->density));

	double leftFlux1 = 0;
	double leftFlux2 = middlePressure*contactDiscontR*contactDiscontR*forwardShockWaveR/forwardV;
	double leftFlux3 = contactDiscontV*leftFlux2;

	double rightFlux1 = upstreamBins2[0]->density*cube(forwardShockWaveR)*(upstreamBins2[0]->U - forwardV)/forwardV;
	double rightFlux2 = cube(forwardShockWaveR)*(upstreamBins2[0]->density*upstreamBins2[0]->U*(upstreamBins2[0]->U - forwardV) + upstreamBins2[0]->pressure)/forwardV;
	double rightFlux3 = cube(forwardShockWaveR)*(upstreamBins2[0]->getEnergy()*(upstreamBins2[0]->U - forwardV) + upstreamBins2[0]->U*upstreamBins2[0]->pressure)/forwardV;

	TracPen(u1, F1, maxSoundSpeed, deltaXi, leftFlux1, rightFlux1);
	TracPen(u2, F2, maxSoundSpeed, deltaXi, leftFlux2, rightFlux2);
	for(int i = 0; i < rgridNumber; ++i){
		u2[i] += tau*downstreamBins2[i]->r*forwardShockWaveR*2*downstreamBins2[i]->pressure*deltaF/forwardV;
	}
	TracPen(u3, F3, maxSoundSpeed, deltaXi, leftFlux3, rightFlux3);

	for(int i = 0; i < rgridNumber; ++i){
		downstreamBins2[i]->u1 = u1[i];
		downstreamBins2[i]->u2 = u2[i];
		downstreamBins2[i]->u3 = u3[i];
	}

	delete[] u1;
	delete[] u2;
	delete[] u3;
	delete[] F1;
	delete[] F2;
	delete[] F3;

}

void Simulation::TracPen(double* u, double* flux, double cs, double deltaXi, double leftFlux, double rightFlux){

	/*u[0] -= tau*((flux[0] + flux[1])/2 - leftFlux)/deltaXi;
	for(int i = 1; i < rgridNumber - 1; ++i){
		u[i] -= tau*0.5*(flux[i+1] - flux[i-1])/deltaXi;
	}
	u[rgridNumber - 1] -= tau*(rightFlux - (flux[rgridNumber - 1] + flux[rgridNumber - 2])/2)/deltaXi;*/
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

	u[0] -= tau*(0.5*(fplus[1] + fminus[1]) - leftFlux)/deltaXi;
	for(int i = 1; i <= rgridNumber-2; ++i){
		u[i] -= tau*0.5*(fplus[i+1] + fminus[i+1] - (fplus[i] + fminus[i]))/deltaXi;
	}
	u[rgridNumber - 1] -= tau*(rightFlux - 0.5*(fplus[rgridNumber - 1] + fminus[rgridNumber - 1]))/deltaXi;

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

void Simulation::moveShockWaves(){
	double oldForwardV = forwardV;
	double oldReverseV = reverseV;
	double oldContactDiscontV = contactDiscontV;

	oldForwardShockWaveR = forwardShockWaveR;
	oldReverseShockWaveR = reverseShockWaveR;
	oldContactDiscontR = contactDiscontR;

	forwardV = upstreamBins2[0]->U + sqrt((gamma + 1)*downstreamBins2[rgridNumber - 1]->pressure/2 + (gamma - 1)*upstreamBins2[0]->pressure/2)/sqrt(upstreamBins2[0]->density);
	if(forwardV <= 0){
		printf("forwardV <= 0\n");
	}
	alertNaNOrInfinity(forwardV,"forwardV");
	reverseV = upstreamBins1[rgridNumber - 1]->U - sqrt((gamma + 1)*downstreamBins1[0]->pressure/2 + (gamma - 1)*upstreamBins1[rgridNumber - 1]->pressure)/sqrt(upstreamBins1[rgridNumber-1]->density);
	if(reverseV <= 0){
	    printf("reverseV <= 0\n");
	}
	alertNaNOrInfinity(reverseV,"reverseV");
	SpaceBin* leftBin = downstreamBins1[rgridNumber - 1];
	SpaceBin* rightBin = downstreamBins2[0];
	contactDiscontV = (leftBin->U*sqrt(leftBin->density) + rightBin->U*sqrt(rightBin->density) + (leftBin->pressure - rightBin->pressure)/(sqrt(gamma*(leftBin->pressure + rightBin->pressure)/2)))/(sqrt(leftBin->density) + sqrt(rightBin->density));
	if(contactDiscontV <= 0){
	    printf("contactDiscontV <= 0\n");
	}
	alertNaNOrInfinity(contactDiscontV, "contactDiscontV");

	forwardShockWaveR *= exp(tau);
	reverseShockWaveR += oldForwardShockWaveR*tau*(reverseV + oldReverseV)*exp(0.5*tau)/(forwardV + oldForwardV);
	contactDiscontR += oldForwardShockWaveR*tau*(contactDiscontV + oldContactDiscontV)*exp(0.5*tau)/(forwardV + oldForwardV);
	time += 2*oldForwardShockWaveR*tau*exp(0.5*tau)/(forwardV + oldForwardV);

	for(int i = 0; i < rgridNumber; ++i){
		upstreamBins2[i]->r = upstreamBins2[i]->xi*forwardShockWaveR;
		downstreamBins2[i]->r = contactDiscontR + downstreamBins2[i]->xi*(forwardShockWaveR - contactDiscontR);
		downstreamBins1[i]->r = contactDiscontR + downstreamBins1[i]->xi*(contactDiscontR - reverseShockWaveR);
		upstreamBins1[i]->r = upstreamBins1[i]->xi*reverseShockWaveR;
	}
	updateBinsVolume();
}

void Simulation::updateBinsVolume(){
    double deltaF = (forwardShockWaveR - contactDiscontR)/rgridNumber;
    double deltaB = (contactDiscontR - reverseShockWaveR)/rgridNumber;

	upstreamBins1[0]->volume = 4*pi*(cube(upstreamBins1[1]->r/2))/3;
	upstreamBins2[0]->volume = 4*pi*(cube(upstreamBins2[1]->r + upstreamBins2[0]->r)/8 - cube(upstreamBins2[0]->r))/3;

	upstreamBins1[rgridNumber - 1]->volume = 4*pi*(cube(upstreamBins1[rgridNumber - 1]->r) - cube(upstreamBins1[rgridNumber - 1]->r + upstreamBins1[rgridNumber - 2]->r)/8)/3;
	upstreamBins2[rgridNumber - 1]->volume = 4*pi*(cube(upstreamBins2[rgridNumber - 1]->r) - cube(upstreamBins2[rgridNumber - 1]->r + upstreamBins2[rgridNumber - 2]->r)/8)/3;

	for(int i = 1; i < rgridNumber - 1; ++i){
		upstreamBins1[i]->volume = pi*(cube(upstreamBins1[i+1]->r + upstreamBins1[i]->r) - cube(upstreamBins1[i]->r + upstreamBins1[i-1]->r))/6;
		upstreamBins2[i]->volume = pi*(cube(upstreamBins2[i+1]->r + upstreamBins2[i]->r) - cube(upstreamBins2[i]->r + upstreamBins2[i-1]->r))/6;
	}

	for(int i = 0; i < rgridNumber; ++i){
		downstreamBins1[i]->volume = 4*pi*(cube(downstreamBins1[i]->r + deltaB/2) - cube(downstreamBins1[i]->r - deltaB/2))/3;
		downstreamBins2[i]->volume = 4*pi*(cube(downstreamBins2[i]->r + deltaF/2) - cube(downstreamBins2[i]->r - deltaF/2))/3;
	}
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

	//tau = 0.00000000000000005*min(deltaF,deltaB)/(rgridNumber*maxSoundSpeed);
	tau = 0.05*1.0/((rgridNumber-1)*maxSoundSpeed);
}

void Simulation::updateParameters(){
	mass = 0;
	downstreamMass1 = 0;
	downstreamMass2 = 0;
	upstreamMass1 = 0;
	upstreamMass2 = 0;
	for(int i = 0; i < rgridNumber; ++i){
		downstreamBins1[i]->temperature = downstreamBins1[i]->pressure*massProton/(downstreamBins1[i]->density*kBoltzman);
		downstreamBins2[i]->temperature = downstreamBins2[i]->pressure*massProton/(downstreamBins2[i]->density*kBoltzman);
		upstreamBins1[i]->temperature = upstreamBins1[i]->pressure*massProton/(upstreamBins1[i]->density*kBoltzman);
		upstreamBins2[i]->temperature = upstreamBins2[i]->pressure*massProton/(upstreamBins2[i]->density*kBoltzman);

		downstreamMass1 += downstreamBins1[i]->density*downstreamBins1[i]->volume;
		downstreamMass2 += downstreamBins2[i]->density*downstreamBins2[i]->volume;
		upstreamMass1 += upstreamBins1[i]->density*upstreamBins1[i]->volume;
		upstreamMass1 += upstreamBins1[i]->density*upstreamBins1[i]->volume;
	}
	mass = downstreamMass1 + downstreamMass2 + upstreamMass1 + upstreamMass2;
}

void Simulation::updateDownstreamValues(){
	double deltaF = forwardShockWaveR - contactDiscontR;
	double deltaB = contactDiscontR - reverseShockWaveR;
	for(int i = 0; i < rgridNumber; ++i){
		downstreamBins1[i]->density = downstreamBins1[i]->u1/(downstreamBins1[i]->r*downstreamBins1[i]->r*deltaB);
		downstreamBins1[i]->U = downstreamBins1[i]->u2/(downstreamBins1[i]->density*downstreamBins1[i]->r*downstreamBins1[i]->r*deltaB);;
		downstreamBins1[i]->pressure = (gamma - 1)*(downstreamBins1[i]->u3/(downstreamBins1[i]->r*downstreamBins1[i]->r*deltaB) - downstreamBins1[i]->density*downstreamBins1[i]->U*downstreamBins1[i]->U/2);

		downstreamBins2[i]->density = downstreamBins2[i]->u1/(downstreamBins2[i]->r*downstreamBins2[i]->r*deltaF);
		downstreamBins2[i]->U = downstreamBins2[i]->u2/(downstreamBins2[i]->density*downstreamBins2[i]->r*downstreamBins2[i]->r*deltaF);;
		downstreamBins2[i]->pressure = (gamma - 1)*(downstreamBins2[i]->u3/(downstreamBins2[i]->r*downstreamBins2[i]->r*deltaF) - downstreamBins2[i]->density*downstreamBins2[i]->U*downstreamBins2[i]->U/2);
	}
}
