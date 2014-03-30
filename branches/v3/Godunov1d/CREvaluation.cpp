#include <time.h>
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

double Simulation::diffussionCoef(int i, int j){
	double p = pgrid[j];
	double B = B0;
	return p*speed_of_light*speed_of_light/(electron_charge*B);
}

//инжекционный член
double Simulation::injection(){
	double pf = pgrid[injectionMomentum];
	double dp = (pgrid[injectionMomentum + 1] - pgrid[injectionMomentum - 1])/2;
	return middleDensity[shockWavePoint]*abs(middleVelocity[shockWavePoint])/(pf*pf*massProton*dp);
	//return 1E-30;
}

//расчет космических лучей

void Simulation::evaluateCR(){

	printf("solve CR\n");
	/*if(shockWavePoint > 0 && shockWavePoint < rgridNumber){
		distributionFunction[shockWavePoint][injectionMomentum] += injection()*deltaT;
	}*/
	double* upper = new double[rgridNumber];
	double* middle = new double[rgridNumber];
	double* lower = new double[rgridNumber];

	double* f = new double[rgridNumber];
 	double* x = new double[rgridNumber];

	double* alpha = new double[rgridNumber];
	double* beta = new double[rgridNumber];

	double* volumeDerivative = new double[rgridNumber];
	double* dtDivdr = new double[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		if(i == rgridNumber - 1){
			volumeDerivative[i] = (middleVelocity[i] - middleVelocity[i-1]);
		} else {
			volumeDerivative[i] = (middleVelocity[i+1] - middleVelocity[i]);
		}
		dtDivdr[i] = deltaT/deltaR[i];
	}

	double deltaLogP = logPgrid[1] - logPgrid[0];
	
	for(int j = 0; j < pgridNumber; ++j){
		double p = pgrid[j];
		for(int i = 0; i < rgridNumber; ++i){
			f[i] = 0;

			double diffCoefLeft = diffussionCoef(i,j);
			double diffCoefRight = diffussionCoef(i+1,j);

			if(i == rgridNumber -1){
				//upper[i] = 0.5*diffCoefRight*dtDivdr[i]/middleDeltaR[i];
				upper[i] = diffCoefRight*dtDivdr[i]/middleDeltaR[i];
			} else {
				//upper[i] =  0.5*diffCoefRight*dtDivdr[i]/middleDeltaR[i+1];
				upper[i] =  diffCoefRight*dtDivdr[i]/middleDeltaR[i+1];
			}
			if(i == 0){
				//lower[i] =  0.5*diffCoefLeft*dtDivdr[i]/middleDeltaR[i];
				lower[i] =  diffCoefLeft*dtDivdr[i]/middleDeltaR[i];
			} else {
				//lower[i] =  0.5*diffCoefLeft*dtDivdr[i]/middleDeltaR[i];
				lower[i] =  diffCoefLeft*dtDivdr[i]/middleDeltaR[i];
			}
			if(i == 0){
				middle[i] = -(upper[i] + lower[i] + 1);
			} else {
				middle[i] = -(upper[i] + lower[i] + 1);
			}

			f[i] = - distributionFunction[i][j];
			/*if(i == 0){
				f[i] = f[i] - upper[i]*distributionFunction[i+1][j] + (upper[i] + lower[i])*distributionFunction[i][j];
			} else if (i == rgridNumber - 1){
				f[i] = f[i] - lower[i]*distributionFunction[i-1][j] + (upper[i] + lower[i])*distributionFunction[i][j];
			} else {
				f[i] = f[i] - lower[i]*distributionFunction[i-1][j] - upper[i]*distributionFunction[i+1][j] + (upper[i] + lower[i])*distributionFunction[i][j];
			}*/

			double derivative = 0;

			if(j == 0){
				derivative = 0;
			} else {
				derivative = (distributionFunction[i][j] - distributionFunction[i][j-1])/deltaLogP;
			}
			f[i] += - volumeDerivative[i]*derivative*dtDivdr[i]/3;

			if((i > 1) && (i < rgridNumber - 1)){
				f[i] += (distributionFunction[i][j]*middleVelocity[i] - distributionFunction[i-1][j]*middleVelocity[i-1])*dtDivdr[i];
 			}
			alertNaNOrInfinity(f[i],"f = NaN");

			if(i == shockWavePoint - 1 && (j == injectionMomentum)){
				f[i] -= injection()*dtDivdr[shockWavePoint-1]*pgrid[j]*pgrid[j]*pgrid[j];
				double dp = (pgrid[j + 1] - pgrid[j - 1])/2;
				injectedParticles += injection()*volume(shockWavePoint - 1)*dtDivdr[shockWavePoint-1]*pgrid[j]*pgrid[j]*dp;
			}
		}

		//не явная
		solveThreeDiagonal(middle, upper, lower, f, x, alpha, beta);
		//явная
		/*for(int i = 0; i < rgridNumber; ++i){
			if(i == 0){
				x[i] = - f[i] + upper[i]*distributionFunction[i+1][j] - (lower[i] + upper[i])*distributionFunction[i][j];
			} else if(i == rgridNumber - 1) {
				x[i] = - f[i] + lower[i]*distributionFunction[i-1][j] - (lower[i] + upper[i])*distributionFunction[i][j];
			} else {
				x[i] = - f[i] + lower[i]*distributionFunction[i-1][j] + upper[i]*distributionFunction[i+1][j] - (lower[i] + upper[i])*distributionFunction[i][j];
			}
		}*/
		
		for(int i = 0; i < rgridNumber; ++i){
			tempDistributionFunction[i][j] = x[i];
			if(x[i] < 0){
				if(abs(x[i]) > 1E-300){
					tempDistributionFunction[i][j] = 0;
					printf("tempDistribution < 0\n");
				} else {
					tempDistributionFunction[i][j] = 0;
				}
			}
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = tempDistributionFunction[i][j];
		}
	}

	evaluateCosmicRayPressure();

	delete[] upper;
	delete[] middle;
	delete[] lower;
	delete[] f;
	delete[] x;
	delete[] alpha;
	delete[] beta;
	delete[] volumeDerivative;
	delete[] dtDivdr;
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
		double momentum = pgrid[j];
		partPressure[j] = momentum*(momentum/sqrt(massProton*massProton + momentum*momentum/(speed_of_light*speed_of_light)))*deltaLogP;
	}
	for(int i = 0; i < rgridNumber; ++i){
		double pressure = 0;
		for(int j = 0; j < pgridNumber; ++j){
			pressure += distributionFunction[i][j]*partPressure[j];
		}
		cosmicRayPressure[i] = 4*pi*pressure;
	}

	delete[] partPressure;
}

void Simulation::changeDistrFunction(){
	for(int j = 0; j < pgridNumber; ++j){
		double p3 = cube(pgrid[j]);
		for(int i = 0; i < rgridNumber; ++i){
			distributionFunction[i][j] *= p3;
		}
	}
}