#include <time.h>
#include "omp.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"
#include "progon.h"

void Simulation::updateDiffusionCoef(){
	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < pgridNumber; ++j){
			double p = pgrid[j];
			//double p = pgrid[0];
			double B = B0;
			/*for(int k = 0; k < kgridNumber; ++k){
				if(kgrid[k] > electron_charge*largeScaleField[i][k]/(speed_of_light*p) ){
					B = largeScaleField[i][k];
					break;
				}
			}*/
			double v = p/sqrt(massProton*massProton + p*p/(speed_of_light*speed_of_light));
			double coef = p*v*speed_of_light/(electron_charge*B);
			double dx = deltaR[i];
			double lambda = coef/speed_of_light;
			if(abs(i - shockWavePoint) < 10 && j >= injectionMomentum){
				if(lambda < dx){
					//printf("lambda < deltaR i = %d j = %d\n",i,j);
				}
			}
			diffusionCoef[i][j] = coef;
		}
	}
}

//инжекционный член
double Simulation::injection(int i){
	double pf = pgrid[injectionMomentum];
	double dp = (pgrid[injectionMomentum + 1] - pgrid[injectionMomentum - 1])/2;
	double xi = pgrid[injectionMomentum]*speed_of_light/(kBoltzman*temperatureIn(i+1));
	double eta = cube(xi)*exp(-xi);
	return (5E-5)*middleDensity[i-1]*abs(middleVelocity[i-1]*middleVelocity[i-1]/speed_of_light)*pf/(massProton*dp*middleDeltaR[i]);
}

//расчет космических лучей

void Simulation::evaluateCR(){

	printf("solve CR\n");
	/*if(shockWavePoint > 0 && shockWavePoint < rgridNumber){
		distributionFunction[shockWavePoint][injectionMomentum] += injection()*deltaT;
	}*/
#pragma omp parallel for
	for(int k = 0; k < pgridNumber; k = k + 1){

	double* upper = new double[rgridNumber+1];
	double* middle = new double[rgridNumber+1];
	double* lower = new double[rgridNumber+1];

	double* f = new double[rgridNumber+1];
 	double* x = new double[rgridNumber+1];

	double* alpha = new double[rgridNumber];
	double* beta = new double[rgridNumber];

		double y = logPgrid[k];
		double p = pgrid[k];
		double gkp = distributionFunction[0][k];
		double gkm=0.0;
		if (k == 0) {
			gkm=0.0;
		} else {
			gkm = tempDistributionFunction[0][k-1];
		}
		if(k == injectionMomentum){
			//gkm = 0;
		}
		//double dx = (rgrid[1] + a)/2;
		double dx = (grid[1] + upstreamR/2)/2;
		double dxp = grid[1] - grid[0];
		double dxm = grid[0] + upstreamR/2;
		double xp = (grid[1] + grid[0])/2;
		double xm = (grid[0] - upstreamR/2)/2;
		middle[0] = (1/(dx))*(diffusionCoef[0][k]/dxp + diffusionCoef[0][k]/dxm) - (1.0/3)*((middleVelocity[1] - middleVelocity[0])/dx)/deltaLogP;
		upper[0] = -(1/(dx))*(diffusionCoef[0][k]/dxp);
		f[0] = - (1.0/3)*((middleVelocity[1] - middleVelocity[0])/dx)*(gkm/deltaLogP);
		for(int i = 1; i < rgridNumber-1; ++i){
			gkp = distributionFunction[i][k];
			if (k == 0) {
				gkm=0;
			} else {
				gkm=tempDistributionFunction[i][k-1];
			}
			//dx = (grid[i+1] - grid[i-1])/2;
			dx = middleDeltaR[i];
			dxp=grid[i+1]-grid[i];
			dxm=grid[i]-grid[i-1];
			xp=(grid[i+1]+grid[i])/2;
			xm=(grid[i]+grid[i-1])/2;
			lower[i-1] = -(1/(dx))*(diffusionCoef[i-1][k]/dxm) - 0.5*middleVelocity[i-1]/dx;
			middle[i] = (1/(dx))*(diffusionCoef[i-1][k]/dxp + diffusionCoef[i][k]/dxm) - (1.0/3)*0.5*((middleVelocity[i+1] - middleVelocity[i-1])/dx)/deltaLogP;
			upper[i] = -(1/(dx))*(diffusionCoef[i][k]/dxp)  + 0.5*middleVelocity[i+1]/dx;

			f[i] = -(1.0/3)*0.5*((middleVelocity[i+1] - middleVelocity[i-1])/dx)*gkm/deltaLogP;

			if(f[i] < 0 ){
				//printf("f[i] < 0\n");
				f[i] = 0;
			}
			if(abs(i-shockWavePoint) < 1 && abs(k - injectionMomentum) < 1){
				double inj = injection(i);
				f[i] += inj;
				//tempDensity[i] -= inj*massProton*deltaLogP;
				//tempMomentum[i] -= inj*massProton*deltaLogP*middleVelocity[i];
				//tempEnergy[i] -= inj*pgrid[k]*speed_of_light*deltaLogP;
				//if(tempDensity[i] < 0){
					//printf("tempDensity[i] < 0 by CR\n");
				//}
				//if(tempEnergy[i] < 0){
					//printf("tempEnergy[i] < 0 by CR\n");
				//}
				//alertNaNOrInfinity(tempDensity[i], "density = NaN");
				//alertNaNOrInfinity(tempMomentum[i], "momentum = NaN");
				//alertNaNOrInfinity(tempEnergy[i], "energy = NaN");
			}
		}
		gkp = distributionFunction[rgridNumber-1][k];
		if (k==0) {
			gkm=0;
		} else {
			gkm = tempDistributionFunction[rgridNumber-1][k-1];
		}
		dx = (grid[rgridNumber-1] - grid[rgridNumber - 2])/2;
		dxm = (grid[rgridNumber-1] - grid[rgridNumber - 2]);
		xp = grid[rgridNumber-1];
		xm = (grid[rgridNumber-1] + grid[rgridNumber-2])/2;
		lower[rgridNumber-2] = -(1/(dx))*(diffusionCoef[rgridNumber-1][k]/dxp);
		middle[rgridNumber-1] = (1/(dx))*(diffusionCoef[rgridNumber-1][k]/dxm) - (1.0/3)*((middleVelocity[rgridNumber-1] - middleVelocity[rgridNumber-2])/dx)/deltaLogP;
		f[rgridNumber-1] = - (1.0/3)*((middleVelocity[rgridNumber-1] - middleVelocity[rgridNumber-2])/dx)*gkm/deltaLogP;
		progon(lower,middle, upper,rgridNumber-1,f,x, alpha, beta);
		if(k == injectionMomentum){
			//outMatrix(lower,middle, upper,rgridNumber-1,f,x);
			//printf("out\n");
		}
		for(int i = 0; i < rgridNumber; ++i){
			//alertNegative(x[i],"tempDistribution < 0");
			alertNaNOrInfinity(x[i],"tempDistribution = NaN");
			tempDistributionFunction[i][k]= x[i];
			if(x[i] < 0){
				tempDistributionFunction[i][k] = 0;
				if(abs(x[i]) > 1E-5){
					printf("tenpDistribution < 0\n");
				}
			}
			if(k < injectionMomentum){
				if(x[i] > 0){
					printf("aaaa\n");
				}
			}
		}
		tempDistributionFunction[rgridNumber][k] =x[rgridNumber-1];
		

	delete[] upper;
	delete[] middle;
	delete[] lower;
	delete[] f;
	delete[] x;
	delete[] alpha;
	delete[] beta;
	}
}

//решение трёх диагональной матрицы
void Simulation::solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x, double* alpha, double* beta){

	alpha[1] = -upper[0]/middle[0];
	beta[1] = f[0]/middle[0];
	for(int i = 2; i < rgridNumber; ++i){
		double temp = lower[i-1]*alpha[i-1] + middle[i-1];
		alpha[i] = -upper[i-1]/temp;
		beta[i] = (f[i-1] - lower[i-1]*beta[i-1])/temp;
	}

	x[rgridNumber - 1] = (f[rgridNumber-1] - lower[rgridNumber-1]*beta[rgridNumber-1])/(lower[rgridNumber-1]*alpha[rgridNumber-1] + middle[rgridNumber-1]);
	alertNaNOrInfinity(x[rgridNumber-1],"x = NaN");
	//alertNegative(x[rgridNumber-1],"x < 0");

	for(int i = rgridNumber - 2; i >= 0; --i){
		x[i] = alpha[i+1]*x[i+1] + beta[i+1];
		alertNaNOrInfinity(x[i],"x = NaN");
		//alertNegative(x[i],"x < 0");
	}
}


//вычисление давления космических лучей

void Simulation::evaluateCosmicRayPressure(){
	double* partPressure = new double[pgridNumber];
	double deltaLogP = logPgrid[1] - logPgrid[0];
	for(int j = 0; j < pgridNumber; ++j){
		double p = pgrid[j];
		double v = p/sqrt(massProton*massProton + p*p/(speed_of_light*speed_of_light));
		partPressure[j] = p*v*deltaLogP;
	}
	for(int i = 0; i < rgridNumber; ++i){
		double pressure = 0;
		double concentration = 0;
		for(int j = 0; j < pgridNumber; ++j){
			if(distributionFunction[i][j] > 0){
				pressure += distributionFunction[i][j]*partPressure[j];
				concentration += distributionFunction[i][j]*deltaLogP;
			}
		}
		//4pi?
		cosmicRayPressure[i] = pressure;
		cosmicRayConcentration[i] = concentration;
	}

	delete[] partPressure;
}