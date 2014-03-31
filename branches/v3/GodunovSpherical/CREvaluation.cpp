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
	double dp = (pgrid[injectionMomentum + 1] + pgrid[injectionMomentum - 1])/2;
	return middleDensity[shockWavePoint]*abs(middleVelocity[shockWavePoint])/(pf*pf*massProton*dp);
	//return 1;
}


//расчет космических лучей

void Simulation::evaluateCR(){
	printf("solve CR\n");
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
		//volumeDerivative[i] = (gridsquare[i+1]*pointVelocity[i+1] - gridsquare[i]*pointVelocity[i]);
		volumeDerivative[i] = (gridsquare[i+1]*middleVelocity[i+1] - gridsquare[i]*middleVelocity[i]);
		dtDivdr[i] = deltaT/deltaR[i];
	}

	double deltaLogP = logPgrid[1] - logPgrid[0];

	/*if(shockWavePoint > 0 && shockWavePoint < rgridNumber){
		distributionFunction[shockWavePoint][injectionMomentum] += injection()*deltaT;
	}*/
	
	for(int j = 0; j < pgridNumber; ++j){
		// -1 ?
		for(int i = 0; i < rgridNumber; ++i){
			//double p = (pgrid[j] + pgrid[j+1])/2;
			double p = pgrid[j];

			double diffCoefLeft = diffussionCoef(i,j);
			double diffCoefRight = diffussionCoef(i+1,j);

			if(i == rgridNumber -1){
				upper[i] = 0.5*gridsquare[i+1]*diffCoefRight*deltaT/(volumeFactor[i]*middleDeltaR[i]);
			} else {
				upper[i] = 0.5*gridsquare[i+1]*diffCoefRight*deltaT/(volumeFactor[i]*middleDeltaR[i+1]);
			}
			if(i == 0){
				lower[i] = 0.5*gridsquare[i]*diffCoefLeft*deltaT/(volumeFactor[i]*middleDeltaR[i]);
			} else {
				lower[i] = 0.5*gridsquare[i]*diffCoefLeft*deltaT/(volumeFactor[i]*middleDeltaR[i]);
			}
			if(i == 0){
				middle[i] = -(upper[i] + lower[i] + 1);
			} else {
				middle[i] = -(upper[i] + lower[i] + 1);
			}
			/*upper[i] = 0;
			middle[i] = -volumeFactor[i];
			lower[i] = 0;*/
			
			f[i] = - distributionFunction[i][j];
			if(i == 0){
				f[i] = f[i] - upper[i]*distributionFunction[i+1][j] + (upper[i] + lower[i])*distributionFunction[i][j];
			} else if (i == rgridNumber - 1){
				f[i] = f[i] - lower[i]*distributionFunction[i-1][j] + (upper[i] + lower[i])*distributionFunction[i][j];
			} else {
				f[i] = f[i] - lower[i]*distributionFunction[i-1][j] - upper[i]*distributionFunction[i+1][j] + (upper[i] + lower[i])*distributionFunction[i][j];
			}

			double derivative = 0;
			if(volumeDerivative[i] >= 0){
				if(j == pgridNumber - 1){
					derivative = 0;
				} else {
					derivative = (distributionFunction[i][j+1] - distributionFunction[i][j])/(3*deltaLogP);
				}
			} else {
				if(j == 0){
					derivative = 0;
				} else {
					derivative = (distributionFunction[i][j] - distributionFunction[i][j-1])/(3*deltaLogP);
				}
			}
			f[i] += - volumeDerivative[i]*derivative*deltaT/volumeFactor[i];

			/*if(middleVelocity[i] > 0){
				if(i > 0){
					f[i] += middleVelocity[i]*(distributionFunction[i][j] - distributionFunction[i-1][j])*deltaT/(volumeFactor[i]*middleDeltaR[i]);
				}
			} else {
				if(i < rgridNumber - 1){
					f[i] += middleVelocity[i]*(distributionFunction[i+1][j] - distributionFunction[i][j])*deltaT/(volumeFactor[i]*middleDeltaR[i+1]);
				}
			}*/
			if(i > 0){
				//f[i] += middleVelocity[i]*middleGrid[i]*middleGrid[i]*(distributionFunction[i][j] - distributionFunction[i-1][j])*deltaT/(volumeFactor[i]);
				//f[i] += (distributionFunction[i][j]*middleVelocity[i] - distributionFunction[i-1][j]*middleVelocity[i-1])*dtDivdr[i];
			}
			alertNaNOrInfinity(f[i],"f = NaN");
			/*if(i == shockWavePoint && j == injectionMomentum){
				f[i] -= injection()*dtDivdr[i];
				double dp = (pgrid[injectionMomentum + 1] - pgrid[injectionMomentum - 1])/2;
				injectedParticles += injection()*volume(shockWavePoint - 1)*deltaT*middleGrid[i]*middleGrid[i]*pgrid[j]*pgrid[j]*dp/volumeFactor[i];
			}*/
		}
		f[rgridNumber-1] -= upper[rgridNumber-1]*distributionFunction[rgridNumber-1][j];

		solveThreeDiagonal(middle, upper, lower, f, x, alpha, beta);
		
		for(int i = 0; i < rgridNumber; ++i){
			tempDistributionFunction[i][j] = x[i];
		}

		/*if(j == injectionMomentum){
			FILE* file = fopen("output/matrix.dat", "w");
			for(int i = 0; i < rgridNumber; ++i){
				fprintf(file,"%d %lf %g %lf %g %g %g %g\n", i, middleGrid[i], lower[i], middle[i], upper[i], f[i], distributionFunction[i][j], tempDistributionFunction[i][j]);
			}
			fclose(file);
		}*/
		if(shockWavePoint > -1){
			//distributionFunction[shockWavePoint-1][injectionMomentum] += injection()*dtDivdr[shockWavePoint-1]*cube(pgrid[injectionMomentum]);
			distributionFunction[shockWavePoint-1][injectionMomentum] += injection()*dtDivdr[shockWavePoint-1];
			double dp = (pgrid[injectionMomentum + 1] - pgrid[injectionMomentum - 1])/2;
			injectedParticles += injection()*volume(shockWavePoint - 1)*dtDivdr[shockWavePoint-1]*pgrid[injectionMomentum]*pgrid[injectionMomentum]*dp;
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
	//alertNaNOrInfinity(x[rgridNumber-1],"x = NaN");
	//alertNegative(x[rgridNumber-1],"x < 0");

	for(int i = rgridNumber - 2; i >= 0; --i){
		x[i] = alpha[i+1]*x[i+1] + beta[i+1];
		//alertNaNOrInfinity(x[i],"x = NaN");
		//alertNegative(x[i],"x < 0");
	}
}


//вычисление давления космических лучей

void Simulation::evaluateCosmicRayPressure(){
	double* partPressure = new double[pgridNumber];
	double deltaLogP = logPgrid[1] - logPgrid[0];
	for(int j = 0; j < pgridNumber; ++j){
		double momentum = pgrid[j];
		partPressure[j] = momentum*momentum*momentum*momentum*(momentum/sqrt(massProton*massProton + momentum*momentum/(speed_of_light*speed_of_light)))*deltaLogP;
	}
	for(int i = 0; i < rgridNumber; ++i){
		double pressure = 0;
		for(int j = 0; j < pgridNumber; ++j){
			//pressure += distributionFunction[i][j]*momentum*momentum*momentum*(momentum/sqrt(massProton*massProton + momentum*momentum/(speed_of_light*speed_of_light)))*(pgrid[j+1] - pgrid[j]);
			pressure += distributionFunction[i][j]*partPressure[j];
		}
		pressure *= 4*pi;
		cosmicRayPressure[i] = pressure;
	}

	delete[] partPressure;
}