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

	double** particleMomentumFlux = new double*[rgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		particleMomentumFlux[i] = new double[pgridNumber + 1];
	}

	double* volumeDerivative = new double[rgridNumber];
	double* dtDivdr = new double[rgridNumber];

	if(currentIteration == changeDistributionParameter) {
		changeDistrFunction();
	}

	for(int i = 0; i < rgridNumber; ++i){
		//volumeDerivative[i] = (gridsquare[i+1]*pointVelocity[i+1] - gridsquare[i]*pointVelocity[i]);
		if(i == rgridNumber - 1){
			volumeDerivative[i] = (middleVelocity[i] - middleVelocity[i-1]);
		} else {
			volumeDerivative[i] = (middleVelocity[i+1] - middleVelocity[i]);
		}
		dtDivdr[i] = deltaT/deltaR[i];
	}

	double deltaLogP = logPgrid[1] - logPgrid[0];

	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j <= pgridNumber; ++j){
			if(j == 0){
				particleMomentumFlux[i][j] = 0;
				continue;
			}
			if(j == pgridNumber){
				particleMomentumFlux[i][j] = 0;
				continue;
			}
			if(j == injectionMomentum){
				if(volumeDerivative[i] > 0){
					particleMomentumFlux[i][j] = volumeDerivative[i]*distributionFunction[i][j-1]/3;
				} else {
					particleMomentumFlux[i][j] = volumeDerivative[i]*distributionFunction[i][j-1]/3;
				}
			} else {
				if(volumeDerivative[i] > 0){
					particleMomentumFlux[i][j] = volumeDerivative[i]*distributionFunction[i][j-1]/3;
				} else {
					particleMomentumFlux[i][j] = volumeDerivative[i]*distributionFunction[i][j-1]/3;
				}
			}
		}
	}
	
	for(int j = 0; j < pgridNumber; ++j){
		double p = pgrid[j];
		//double p3 = p*p*p;
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
			/*if(volumeDerivative[i] > 0){
				if((j >= pgridNumber - 2) || (j == 0)){
					derivative = 0;
					distrFunDerivative[j] = 0;
					distrFunDerivative2[j] = 0;
				} else {
					//derivative = (distributionFunction[i][j+1] - distributionFunction[i][j])/(3*deltaLogP);
				//}
					double a,b,c,d, p1,p2,p3,p4, y1,y2,y3,y4;
					p1 = pgrid[j-1];
					p2 = pgrid[j];
					p3 = pgrid[j+1];
					p4 = pgrid[j+2];
					y1 = distributionFunction[i][j-1];
					y2 = distributionFunction[i][j];
					y3 = distributionFunction[i][j+1];
					y4 = distributionFunction[i][j+2];
	
					//a = y1/((p1-p4)*(p1-p2)*(p1-p3)) + y2/((p2-p1)*(p2-p3)*(p2-p4)) + y3/((p3-p1)*(p3-p2)*(p3-p4)) + y4/((p4-p1)*(p4-p2)*(p4-p3));
					//b = y1*(-p2-p3-p4)/((p1-p4)*(p1-p2)*(p1-p3)) + y2*(-p1-p3-p4)/((p2-p1)*(p2-p3)*(p2-p4)) + y3*(-p1-p2-p4)/((p3-p1)*(p3-p2)*(p3-p4)) + y4*(-p1-p2-p3)/((p4-p1)*(p4-p2)*(p4-p3));
					//c = y1*(p2*p3+p2*p4+p3*p4)/((p1-p4)*(p1-p2)*(p1-p3)) + y2*(p1*p3 + p1*p4 + p3*p4)/((p2-p1)*(p2-p3)*(p2-p4)) + y3*(p1*p2 + p1*p4 + p2*p4)/((p3-p1)*(p3-p2)*(p3-p4)) + y4*(p1*p2 + p1*p3 + p2*p3)/((p4-p1)*(p4-p2)*(p4-p3));
					//d = y1*(-p2*p3*p4)/((p1-p4)*(p1-p2)*(p1-p3)) + y2*(-p1*p3*p4)/((p2-p1)*(p2-p3)*(p2-p4)) + y3*(-p1*p2*p4)/((p3-p1)*(p3-p2)*(p3-p4)) + y4*(-p1*p2*p3)/((p4-p1)*(p4-p2)*(p4-p3));

					//derivative = 3*a*p2*p2 + 2*b*p2 + c;

					//derivative = (-3*y2 + 4*y3 - y4)/(2*deltaLogP);


					derivative = (y3 - y2)/deltaLogP;

					if(i == shockWavePoint){
						distrFunDerivative[j] = derivative;
						distrFunDerivative2[j] = (y3 - y2)/deltaLogP;
					}

					if(j == injectionMomentum){
						derivative = (y3 - y1)/(2*deltaLogP);
					}

					//f[i] += volumeDerivative[i]*y3*dtDivdr[i]/(3*deltaLogP);
					//middle[i] += volumeDerivative[i]*dtDivdr[i]/(3*deltaLogP);
				}
			} else {
				if((j < 2) || (j == rgridNumber - 1)){
					derivative = 0;
					distrFunDerivative[j] = 0;
					distrFunDerivative2[j] = 0;
				} else {
					double a,b,c,d, p1,p2,p3,p4, y1,y2,y3,y4;
					p1 = pgrid[j-2];
					p2 = pgrid[j-1];
					p3 = pgrid[j];
					p4 = pgrid[j+1];
					y1 = distributionFunction[i][j-2];
					y2 = distributionFunction[i][j-1];
					y3 = distributionFunction[i][j];
					y4 = distributionFunction[i][j-1];
	
					//a = y1/((p1-p4)*(p1-p2)*(p1-p3)) + y2/((p2-p1)*(p2-p3)*(p2-p4)) + y3/((p3-p1)*(p3-p2)*(p3-p4)) + y4/((p4-p1)*(p4-p2)*(p4-p3));
					//b = y1*(-p2-p3-p4)/((p1-p4)*(p1-p2)*(p1-p3)) + y2*(-p1-p3-p4)/((p2-p1)*(p2-p3)*(p2-p4)) + y3*(-p1-p2-p4)/((p3-p1)*(p3-p2)*(p3-p4)) + y4*(-p1-p2-p3)/((p4-p1)*(p4-p2)*(p4-p3));
					//c = y1*(p2*p3+p2*p4+p3*p4)/((p1-p4)*(p1-p2)*(p1-p3)) + y2*(p1*p3 + p1*p4 + p3*p4)/((p2-p1)*(p2-p3)*(p2-p4)) + y3*(p1*p2 + p1*p4 + p2*p4)/((p3-p1)*(p3-p2)*(p3-p4)) + y4*(p1*p2 + p1*p3 + p2*p3)/((p4-p1)*(p4-p2)*(p4-p3));
					//d = y1*(-p2*p3*p4)/((p1-p4)*(p1-p2)*(p1-p3)) + y2*(-p1*p3*p4)/((p2-p1)*(p2-p3)*(p2-p4)) + y3*(-p1*p2*p4)/((p3-p1)*(p3-p2)*(p3-p4)) + y4*(-p1*p2*p3)/((p4-p1)*(p4-p2)*(p4-p3));

					//derivative = 3*a*p3*p3 + 2*b*p3 + c;

					//derivative = -(-3*y3 + 4*y2 - y1)/(2*deltaLogP);

					derivative = (y3 - y2)/deltaLogP;

					if(i == shockWavePoint){
						distrFunDerivative[j] = derivative;
						distrFunDerivative2[j] = (y3 - y2)/deltaLogP;
					}

					if(j == injectionMomentum){
						derivative = (y4 - y2)/(2*deltaLogP);
					}

					//f[i] += volumeDerivative[i]*y2*dtDivdr[i]/(3*deltaLogP);
					//middle[i] += volumeDerivative[i]*dtDivdr[i]/(3*deltaLogP);
				}
			}*/

			if(j == 0){
				derivative = 0;
			} else {
				derivative = (distributionFunction[i][j] - distributionFunction[i][j-1])/deltaLogP;
			}
			f[i] += - volumeDerivative[i]*derivative*dtDivdr[i]/3;
			//f[i] += - volumeDerivative[i]*derivative*dtDivdr[i]*p/3;

			//f[i] -= (particleMomentumFlux[i][j+1] - particleMomentumFlux[i][j])*dtDivdr[i]/deltaLogP;

			if((i > 1) && (i < rgridNumber - 1)){
				//middle[i] -= middleVelocity[i]*dtDivdr[i];
				//lower[i] += middleVelocity[i]*dtDivdr[i];
				//f[i] += middleVelocity[i]*(distributionFunction[i][j] - distributionFunction[i-1][j])*dtDivdr[i];
				//f[i] += middleVelocity[i]*(3*distributionFunction[i][j] - 4*distributionFunction[i-1][j] + 4*distributionFunction[i-2][j])*0.5*dtDivdr[i];
				//f[i] += (distributionFunction[i][j]*middleVelocity[i] - distributionFunction[i-1][j]*middleVelocity[i-1])*dtDivdr[i];

				if(currentIteration < changeDistributionParameter){
					f[i] += 0.5*(middleVelocity[i]+middleVelocity[i-1])*(distributionFunction[i][j] - distributionFunction[i-1][j])*dtDivdr[i];
				} else {
					f[i] += (distributionFunction[i][j]*middleVelocity[i] - distributionFunction[i-1][j]*middleVelocity[i-1])*dtDivdr[i];
				}
 			}
			alertNaNOrInfinity(f[i],"f = NaN");
			//alertNegative(-f[i], "f > 0");

			if(i == shockWavePoint - 1 && (j == injectionMomentum) && (currentIteration < 5000000000)){
				//f[i] -= injection()*dtDivdr[shockWavePoint-1]*pgrid[j]*pgrid[j]*pgrid[j];
				//f[i] -= injection()*dtDivdr[shockWavePoint-1];
				if(currentIteration < changeDistributionParameter){
					f[i] -= injection()*dtDivdr[shockWavePoint-1];
				} else {
					f[i] -= injection()*dtDivdr[shockWavePoint-1]*pgrid[j]*pgrid[j]*pgrid[j];
				}
				double dp = (pgrid[j + 1] - pgrid[j - 1])/2;
				injectedParticles += injection()*volume(shockWavePoint - 1)*dtDivdr[shockWavePoint-1]*pgrid[j]*pgrid[j]*dp;
			}
		}
		//f[rgridNumber-1] -= upper[rgridNumber-1]*distributionFunction[rgridNumber-1][j];

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

		/*if(j == injectionMomentum + 1){
			FILE* file = fopen("output/matrix.dat", "w");
			for(int i = 0; i < rgridNumber; ++i){
				fprintf(file,"%d %lf %g %lf %g %g %g %g\n", i, middleGrid[i], lower[i], middle[i], upper[i], f[i], distributionFunction[i][j], tempDistributionFunction[i][j]);
			}
			fclose(file);
		}*/
	}

	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = tempDistributionFunction[i][j];
			/*if(distributionFunction[i][j] < 0){
				printf("distributionFunction < 0\n");
			}*/
		}
	}

	/*if(shockWavePoint > -1){
		//distributionFunction[shockWavePoint-1][injectionMomentum] += injection()*dtDivdr[shockWavePoint-1]*cube(pgrid[injectionMomentum]);
		distributionFunction[shockWavePoint-1][injectionMomentum] += injection()*dtDivdr[shockWavePoint-1];
		double dp = (pgrid[injectionMomentum + 1] - pgrid[injectionMomentum - 1])/2;
		injectedParticles += injection()*volume(shockWavePoint - 1)*dtDivdr[shockWavePoint-1]*pgrid[injectionMomentum]*pgrid[injectionMomentum]*dp;
	}*/

	evaluateCosmicRayPressure();

	for(int i = 0; i < rgridNumber; ++i){
		delete[] particleMomentumFlux[i];
	}
	delete[] particleMomentumFlux;

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
		partPressure[j] = momentum*momentum*momentum*momentum*(momentum/sqrt(massProton*massProton + momentum*momentum/(speed_of_light*speed_of_light)))*deltaLogP;
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