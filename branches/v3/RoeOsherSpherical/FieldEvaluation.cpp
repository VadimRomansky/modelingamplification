#include "math.h"
#include "stdio.h"
#include <stdlib.h>
#include <omp.h>
#include "simulation.h"
#include "constants.h"
#include "util.h"
#include "complex.h"

void Simulation::evaluateField(){
	//printf("evaluating magnetic field\n");

	//growthRate();

	double fullMaxRate = 0;
	int i;
	#pragma omp parallel for private(i)
	for(i = 1; i < rgridNumber; ++i){
		int maxRateK = 0;
		maxRate[i] = 0;
		for(int k = 0; k < kgridNumber; ++k){
			if(middleVelocity[i] > 0 && middleVelocity[i-1] > 0){
				tempMagneticField[i][k] = magneticField[i][k];
				double Ualpha = power(middleVelocity[i],1.5);
				double prevUalpha = power(middleVelocity[i-1],1.5);
				double z = magneticField[i][k]*Ualpha*middleGrid[i]*middleGrid[i];
				double prevZ = magneticField[i-1][k]*prevUalpha*middleGrid[i-1]*middleGrid[i-1];
				double nextZ = z;
				if(i < rgridNumber-1) {
					nextZ = magneticField[i+1][k]*power(middleVelocity[i+1],1.5)*middleGrid[i+1]*middleGrid[i+1];
				}
				tempMagneticField[i][k] -= deltaT*(z - prevZ)/(sqrt(middleVelocity[i])*sqr(middleGrid[i])*deltaR[i]);
			} else if(middleVelocity[i] < 0 && middleVelocity[i-1] < 0){
				double Ualpha = power(-middleVelocity[i],1.5);
				double z = magneticField[i][k]*Ualpha*middleGrid[i]*middleGrid[i];
				double prevZ = magneticField[i-1][k]*power(-middleVelocity[i-1],1.5)*middleGrid[i-1]*middleGrid[i-1];
				tempMagneticField[i][k] += deltaT*(z - prevZ)/(sqrt(-middleVelocity[i])*sqr(middleGrid[i])*deltaR[i]);
			} else {
				tempMagneticField[i][k] += -deltaT*1.5*((sqr(middleGrid[i])*middleVelocity[i]*magneticField[i][k] - sqr(middleGrid[i-1])*middleVelocity[i-1]*magneticField[i-1][k])/(sqr(middleGrid[i])*deltaR[i])) + deltaT*0.5*0.5*(middleVelocity[i]+middleVelocity[i-1])*(magneticField[i][k] - magneticField[i-1][k])/deltaR[i];
			}

			tempMagneticField[i][k] +=  deltaT*growth_rate[i][k]*magneticField[i][k];
			if(growth_rate[i][k] > maxRate[i]){
				maxRateK = k;
				maxRate[i] = growth_rate[i][k];
			}
			//tempMagneticField[i][k] = magneticField[i][k] + deltaT*(- 1.5*(magneticField[i][k]*middleVelocity[i] - magneticField[i-1][k]*middleVelocity[i-1])/middleDeltaR[i] + 0.5*(delta*middleVelocity[i]+(1-delta)*middleVelocity[i-1])*(magneticField[i][k] - magneticField[i-1][k])/middleDeltaR[i] + growth_rate[i][k]);
			alertNaNOrInfinity(tempMagneticField[i][k], "magnetic field = NaN");
			if(tempMagneticField[i][k] < 0){
				printf("magneticField < 0\n");
				tempMagneticField[i][k] = 0;
			}
		}
		if(maxRate[i] > fullMaxRate){
			fullMaxRate = maxRate[i];
		}
	}

}

void Simulation::evaluateCRFlux(){
	int i;
	#pragma omp for
		for(i = 0; i < rgridNumber; ++i){
			for(int j = 1; j < pgridNumber; ++j){
				crflux[i][j] = - electron_charge*(diffusionCoef[i][j]*(distributionFunction[i+1][j] - distributionFunction[i][j])/deltaR[i])*deltaLogP;
			}
		}
}

void Simulation::growthRate(){
	int i;
    #pragma omp parallel for private(i)
		for(i = 0; i < rgridNumber; ++i){
			if(i < shockWavePoint){
				for(int k = 0; k < kgridNumber; ++k){
					growth_rate[i][k] = 0;
				}
				continue;
			}
			double J = 0;
			for(int j = 0; j < pgridNumber; ++j){
				J += crflux[i][j];
			}
			integratedFlux[i] = J;

			if(J == 0){
				for(int k = 0; k < kgridNumber; ++k){
					growth_rate[i][k] = 0;
				}
				continue;
			}

			double Bls = B0;
			for(int k = 0; k < kgridNumber; ++k){
				Bls = sqrt(4*pi*magneticField[i][k]*kgrid[k]*deltaLogK + Bls*Bls);
				alertNaNOrInfinity(Bls, "Bls = NaN");
                double kc = abs2(4*pi*J/(speed_of_light*Bls));
				double Va = Bls/sqrt(4*pi*middleDensity[i]);
				Complex A1;
				Complex A2;
				for(int j = goodMomentum; j < pgridNumber; ++j){
					double z = kgrid[k]*speed_of_light*pgrid[j]/(electron_charge*Bls);
					alertNaNOrInfinity(z, "z = NaN");
					Complex sigma1;
					Complex sigma2;
                    if( abs2(z - 1) < 0.00001){
						sigma1 = 3/2;
					} else if(z > 1) {
                        sigma1 = (1.5/sqr(z)) + 0.75*(1 - 1/(sqr(z)))*log(abs2((z+1)/(z-1)))/z;
					} else if(0.01 < z) {
                        sigma1 = (1.5/sqr(z)) + 0.75*(1 - 1/(sqr(z)))*log(abs2((z+1)/(z-1)))/z;
					} else {
						sigma1 = 1 + 0.2*sqr(z);
					}

					if(z > 1){
						sigma1 = sigma1 + Complex(0,-(3*pi/(4*z))*(1 - 1/sqr(z)));
					}
					alertNaNOrInfinity(sigma1.real, "sigma = NaN");
					alertNaNOrInfinity(sigma1.im, "sigma = NaN");
					sigma2 = sigma1.conjugate();

					A1 = A1 + sigma1*crflux[i][j];

					alertNaNOrInfinity(A1.real, "A = NaN");
					alertNaNOrInfinity(A1.im, "A = NaN");
					A2 = A2 + sigma2*crflux[i][j];
				}	
				Complex complex1 = Complex(1);
				Complex b1 = (complex1 - (A1/J - 1)*(kc/kgrid[k]))*sqr(kgrid[k]*Va);
				Complex b2 = (complex1 + (A2/J - 1)*(kc/kgrid[k]))*sqr(kgrid[k]*Va);
				//double alpha = 1.5;
				double alpha = 0;
				Complex d1 = Complex(0, -1)*((A1*0.5/J) + 1.5)*(kgrid[k]*kc)*alpha/(4*pi*middleDensity[i]);
				Complex d2 = Complex(0, 1)*((A2*0.5/J) + 1.5)*(kgrid[k]*kc)*alpha/(4*pi*middleDensity[i]);
				Complex G1p = (csqrt(d1*d1 +b1*4) - d1)/2;
				Complex G1m = (csqrt(d1*d1 +b1*4) + d1)/(-2);

				Complex G2p = (csqrt(d2*d2 +b2*4) - d2)/2;
				Complex G2m = (csqrt(d2*d2 +b2*4) + d2)/(-2);

				double rate = max2(max2(G1p.im, G1m.im), max2(G2p.im, G2m.im));
				//alertNaNOrInfinity(rate, "rate = NaN");
				if(rate > 0){
					rate *= 2;
				} else {
					rate = 0;
				}

				if((rate != rate) || (0*rate != 0*rate)){
					printf("rate = NaN\n");
				}

				growth_rate[i][k] = rate;
				alertNaNOrInfinity(growth_rate[i][k], "growth_rate[i][k] = NaN");
			}
		}
	
}
