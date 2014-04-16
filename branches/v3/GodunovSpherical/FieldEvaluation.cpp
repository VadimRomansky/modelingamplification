#include "simulation.h"
#include "constants.h"
#include "util.h"
#include "complex.h"

void Simulation::evauateField(){

}

void Simulation::evaluateCRFlux(double* crflux){
}

void Simulation::growthRate(double** crflux, double** growth_rate){
	for(int i = 0; i < rgridNumber; ++i){
		double J = 0;
		for(int j = 0; j < pgridNumber; ++j){
			J = crflux[i][j]*cube(pgrid[j])*deltaLogP;
		}


		for(int k = 0; k < kgridNumber; ++i){
			double Bls = sqrt(4*pi*magneticField[i][k] + B0*B0);
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

			growth_rate[i][k] = 
		}
	}
}