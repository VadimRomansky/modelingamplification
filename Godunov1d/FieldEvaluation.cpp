#include "simulation.h"
#include "constants.h"
#include "util.h"
#include "complex.h"

void Simulation::evaluateField(){
	evaluateCRFlux();
	growthRate();

	for(int i = 1; i < rgridNumber; ++i){
		for(int k = 0; k < kgridNumber; ++k){
			tempMagneticField[i][k] = magneticField[i][k] + deltaT*(- 1.5*magneticField[i][k]*(middleVelocity[i] - middleVelocity[i-1])/middleDeltaR[i] - middleVelocity[i]*(magneticField[i-1][k] - magneticField[i][k])/middleDeltaR[i] + growth_rate[i][k]);
			if(tempMagneticField[i][k] < 0){
				tempMagneticField[i][k] = 0;
			}
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		for(int k = 0; k < kgridNumber; ++k){
			magneticField[i][k] = tempMagneticField[i][k];
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		double magneticEnergy = 0;
		for(int k = 0; k < rgridNumber; ++k){
			magneticEnergy += magneticField[i][k];
		}
		magneticInductionSum[i] = sqrt(4*pi*magneticEnergy + B0*B0);
	}
}

void Simulation::evaluateCRFlux(){
	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < pgridNumber; ++j){
			crflux[i][j] = electron_charge*diffusionCoef(i,pgrid[j])*(distributionFunction[i+1][j] - distributionFunction[i][j])/deltaR[i];
		}
	}
}

void Simulation::growthRate(){
	for(int i = 0; i < rgridNumber; ++i){
		double J = 0;
		for(int j = 0; j < pgridNumber; ++j){
			J = crflux[i][j]*cube(pgrid[j])*deltaLogP;
		}

		double Bls = B0;
		for(int k = 0; k < kgridNumber; ++k){
			Bls = sqrt(4*pi*magneticField[i][k]*dk + Bls*Bls);
			double kc = abs(4*pi*J/(speed_of_light*Bls));
			double Va = Bls/sqrt(4*pi*middleDensity[i]);
			Complex A1;
			Complex A2;
			for(int j = 0; j < pgridNumber; ++j){
				double z = kgrid[k]*speed_of_light*pgrid[j]/(electron_charge*Bls);
				Complex sigma1;
				Complex sigma2;
				if( z == 1){
					sigma1 = 3/2;
				} else if(z > 1) {
					sigma1 = (1.5/sqr(z)) + 0.75*(1 - 1/(sqr(z)))*log(abs((z+1)/(z-1)))/z;
				} else if(0.01 < z) {
					sigma1 = (1.5/sqr(z)) + 0.75*(1 - 1/(sqr(z)))*log(abs((z+1)/(z-1)))/z;
				} else {
					sigma1 = 1 + 0.2*sqr(z);
				}

				if(z > 1){
					sigma1 = sigma1 + Complex(0,-(3*pi/(4*z))*(1 - 1/sqr(z)));
				}

				sigma2 = sigma1.conjugate();

				A1 = A1 + sigma1*(crflux[i][j]*cube(pgrid[j])*deltaLogP);
				A2 = A2 + sigma2*(crflux[i][j]*cube(pgrid[j])*deltaLogP);
			}	
			Complex b1 = (A1*electron_charge/J - 1)*sqr(kgrid[k]*Va)*(1 - (kc/kgrid[k]));
			Complex b2 = (A2*electron_charge/J - 1)*sqr(kgrid[k]*Va)*(1 + (kc/kgrid[k]));
			double alpha = 1.5;
			Complex d1 = Complex(0, -1)*((A1*0.5*electron_charge/J) + 1.5)*(kgrid[k]*kc)*alpha/(4*pi*middleDensity[i]);
			Complex d2 = Complex(0, 1)*((A2*0.5*electron_charge/J) + 1.5)*(kgrid[k]*kc)*alpha/(4*pi*middleDensity[i]);

			Complex G1p = (csqrt(d1*d1 +b1*4) - d1)/2;
			Complex G1m = (csqrt(d1*d1 +b1*4) + d1)/(-2);

			Complex G2p = (csqrt(d2*d2 +b2*4) - d2)/2;
			Complex G2m = (csqrt(d2*d2 +b2*4) + d2)/(-2);

			double rate = max2(max2(G1p.im, G1m.im), max2(G2p.im, G2m.im));
			if(rate > 0){
				rate *= 2;
			} else {
				rate = 0;
			}

			growth_rate[i][k] = rate*magneticField[i][k];
		}
	}
}