#include "stdafx.h"
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

Simulation::Simulation(){
	cycleBound = true;
	time = 0;
}

Simulation::~Simulation(){
	delete[] grid;
	delete[] pointDensity;
	delete[] pointVelocity;
	delete[] pointPressure;
	delete[] middleDensity;
	delete[] middleVelocity;
	delete[] middlePressure;
}

void Simulation::initializeProfile(){
	downstreamR = 0;
	grid = new double[rgridNumber +1];
	pointDensity = new double[rgridNumber + 1];
	pointVelocity = new double[rgridNumber + 1];
	pointPressure = new double[rgridNumber + 1];
	middleDensity = new double[rgridNumber];
	middleVelocity = new double[rgridNumber];
	middlePressure = new double[rgridNumber];
	deltaR = (upstreamR - downstreamR)/rgridNumber;
	double r = downstreamR;
	double pressure0 = density0*kBoltzman*temperature/massProton;
	for(int i = 0; i < rgridNumber + 1; ++i){
		grid[i] = r;
		r += deltaR;
		switch(simulationType){
		case 1 :
			middleDensity[i] = density0;
			if(i > rgridNumber/10 && i < 2*rgridNumber/10){
				middleVelocity[i] = U0;
			} else {
				middleVelocity[i] = 0;
			}
			middlePressure[i] = pressure0;
			break;
		case 2 :
			middleDensity[i] = density0 + density0*0.01*sin(i*10*2*pi/rgridNumber);
			middleVelocity[i] = sqrt(gamma*pressure0/density0)*density0*0.01*sin(i*10*2*pi/rgridNumber)/density0;
			middlePressure[i] = pressure0 + (gamma*pressure0/density0)*density0*0.01*sin(i*10*2*pi/rgridNumber);
			break;
		default:
			middleDensity[i] = density0;
			middleVelocity[i] = 0;
			if(i <= rgridNumber/100){
				middlePressure[i] = 100000000*pressure0;
			} else {
				middlePressure[i] = pressure0;
			}
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		pointDensity[i] = middleDensity[i];
		pointVelocity[i] = middleVelocity[i];
		pointPressure[i] = middlePressure[i];
	}

	if(cycleBound){
		pointDensity[rgridNumber] = pointDensity[0];
		pointVelocity[rgridNumber] = pointVelocity[0];
		pointPressure[rgridNumber] = pointPressure[0];
	} else {
		pointDensity[rgridNumber] = pointDensity[rgridNumber - 1];
		pointVelocity[rgridNumber] = pointVelocity[rgridNumber - 1];
		pointPressure[rgridNumber] = pointPressure[rgridNumber - 1];
	}
}

void Simulation::simulate(){
	FILE* outFile = fopen("./output/tamc_radial_profile.dat","w");
	FILE* outIteration = fopen("./output/iterations.dat","w");
	fclose(outIteration);
	printf("initialization\n");
	initializeProfile();
	output(outFile,this);
	fclose(outFile);
	updateMaxSoundSpeed();
	for(int i = 0; i < iterationNumber; ++i){
		printf("iteration � %d\n", i);
		printf("time = %lf\n", time);
		printf("solving\n");
		time = time + deltaT;
		evaluateHydrodynamic();
		//updateValues();
		updateMaxSoundSpeed();
		//updateParameters();
		if(i % 10 == 0){
			printf("outputing\n");
			outFile = fopen("./output/tamc_radial_profile.dat","a");
			output(outFile, this);
			fclose(outFile);
			outIteration = fopen("./output/iterations.dat","a");
			fprintf(outIteration, "%d %lf %lf\n", i, time, mass);
			fclose(outIteration);
		}
	}
}

void Simulation::evaluateHydrodynamic(){
	solveDiscontinious();

	for(int i = 0; i < rgridNumber; ++i){
		middleDensity[i] -= deltaT*(densityFlux(i+1) - densityFlux(i))/deltaR;
		double middleMomentum = momentum(i) - deltaT*(momentumFlux(i+1) - momentumFlux(i))/deltaR;
		double middleEnergy = energy(i) - deltaT*(energyFlux(i+1) - energyFlux(i))/deltaR;
		middleVelocity[i] = middleMomentum/middleDensity[i];
		//middleDensity[i] -= deltaT*(densityFlux(i+1) - densityFlux(i))/deltaR;
		middlePressure[i] = (middleEnergy - middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2)*(gamma - 1);

		if(middleDensity[i] <= epsilon*density0){
			middleDensity[i] = epsilon*density0;
			middleVelocity[i] = 0;
			middlePressure[i] = middleDensity[i]*kBoltzman*temperature/massProton;
		}
		if(middlePressure[i] <= 0){
			middlePressure[i] = epsilon*middleDensity[i]*kBoltzman*temperature/massProton;
		}
	}
}

void Simulation::solveDiscontinious(){
	for(int i = 1; i < rgridNumber; ++i){
		double p = pointPressure[i];
		double p1 = middlePressure[i-1];
		double p2 = middlePressure[i];
		double rho1 = middleDensity[i-1];
		double rho2 = middleDensity[i];
		double u1 = middleVelocity[i-1];
		double u2 = middleVelocity[i];
		double u;
		double R1;
		double R2;
		double alpha1;
		double alpha2;
		
		successiveApproximatonPressure(p, u, R1, R2, alpha1, alpha2, p1, p2, u1, u2, rho1, rho2);

		bool isLeftShockWave = (p < p1);
		bool isRightShockWave = (p < p2);

		double D1;
		double D2;

		double c1 = sqrt(gamma*p1/rho1);
		double c2 = sqrt(gamma*p2/rho2);
		if(isLeftShockWave){
			D1 = u1 - alpha1/rho1;
		} else {
			D1 = u1 - c1;
		}

		if(isRightShockWave){
			D2 = u2 + alpha2/rho2;
		} else {
			D2 = u2 + c2;
		}

		if(D1 > 0){
			pointDensity[i] = middleDensity[i - 1];
			pointVelocity[i] = middleVelocity[i - 1];
			pointPressure[i] = middlePressure[i - 1];
		} else if(D2 < 0){
			pointDensity[i] = middleDensity[i];
			pointVelocity[i] = middleVelocity[i];
			pointPressure[i] = middlePressure[i];
		} else {
			if(u > 0){
				if(isLeftShockWave){
					pointDensity[i] = R1;
					pointVelocity[i] = u;
					pointPressure[i] = p;
				} else {
					double D3 = u - c1 - (gamma - 1)*(u1 - u)/2;
					if(D3 < 0){
						pointDensity[i] = R1;
						pointVelocity[i] = u;
						pointPressure[i] = p;
					} else {
						pointVelocity[i] = (gamma - 1)*middleVelocity[i-1]/(gamma + 1) + 2*c1/(gamma + 1);
						pointPressure[i] = middlePressure[i-1]*power(pointVelocity[i]/c1, 2*gamma/(gamma - 1));
						pointDensity[i] = gamma*pointPressure[i]/(pointVelocity[i]*pointVelocity[i]);
					}
				}
			} else {
				if(isRightShockWave){
					pointDensity[i] = R2;
					pointVelocity[i] = u;
					pointPressure[i] = p;
				} else {
					double D3 = u + c2 - (gamma - 1)*(u2 - u)/2;
					if(D3 < 0){
						/////????????
						pointVelocity[i] = (gamma - 1)*middleVelocity[i]/(gamma + 1) + 2*c2/(gamma + 1);
						pointPressure[i] = middlePressure[i]*power(pointVelocity[i]/c2, 2*gamma/(gamma - 1));
						pointDensity[i] = gamma*pointPressure[i]/(pointVelocity[i]*pointVelocity[i]);
					} else {
						pointDensity[i] = R2;
						pointVelocity[i] = u;
						pointPressure[i] = p;
					}
				}
			}
		}
	}

	if(cycleBound){
		double p = pointPressure[0];
		double p1 = middlePressure[rgridNumber-1];
		double p2 = middlePressure[0];
		double rho1 = middleDensity[rgridNumber-1];
		double rho2 = middleDensity[0];
		double u1 = middleVelocity[rgridNumber-1];
		double u2 = middleVelocity[0];
		double u;
		double R1;
		double R2;
		double alpha1;
		double alpha2;
		
		successiveApproximatonPressure(p, u, R1, R2, alpha1, alpha2, p1, p2, u1, u2, rho1, rho2);

		double c1 = sqrt(gamma*p1/rho1);
		double c2 = sqrt(gamma*p2/rho2);

		bool isLeftShockWave = (p < p1);
		bool isRightShockWave = (p < p2);

		double D1;
		double D2;

		if(isLeftShockWave){
			D1 = u1 - alpha1/rho1;
		} else {
			D1 = u1 - c1;
		}

		if(isRightShockWave){
			D2 = u2 + alpha2/rho2;
		} else {
			D2 = u2 + c2;
		}

		if(D1 > 0){
			pointDensity[0] = middleDensity[rgridNumber - 1];
			pointVelocity[0] = middleVelocity[rgridNumber - 1];
			pointPressure[0] = middlePressure[rgridNumber - 1];
		} else if(D2 < 0){
			pointDensity[0] = middleDensity[0];
			pointVelocity[0] = middleVelocity[0];
			pointPressure[0] = middlePressure[0];
		} else {
			if(u > 0){
				if(isLeftShockWave){
					pointDensity[0] = R1;
					pointVelocity[0] = u;
					pointPressure[0] = p;
				} else {
					double D3 = u - c1 - (gamma - 1)*(u1 - u)/2;
					if(D3 < 0){
						pointDensity[0] = R1;
						pointVelocity[0] = u;
						pointPressure[0] = p;
					} else {
						pointVelocity[0] = (gamma - 1)*middleVelocity[rgridNumber - 1]/(gamma + 1) + 2*c1/(gamma + 1);
						pointPressure[0] = middlePressure[rgridNumber - 1]*power(pointVelocity[0]/c1, 2*gamma/(gamma - 1));
						pointDensity[0] = gamma*pointPressure[0]/(pointVelocity[0]*pointVelocity[0]);
					}
				}
			} else {
				if(isRightShockWave){
					pointDensity[0] = R2;
					pointVelocity[0] = u;
					pointPressure[0] = p;
				} else {
					double D3 = u + c2 - (gamma - 1)*(u2 - u)/2;
					if(D3 < 0){
						/////????????
						pointVelocity[0] = (gamma - 1)*middleVelocity[0]/(gamma + 1) + 2*c2/(gamma + 1);
						pointPressure[0] = middlePressure[0]*power(pointVelocity[0]/c2, 2*gamma/(gamma - 1));
						pointDensity[0] = gamma*pointPressure[0]/(pointVelocity[0]*pointVelocity[0]);
					} else {
						pointDensity[0] = R2;
						pointVelocity[0] = u;
						pointPressure[0] = p;
					}
				}
			}
		}

		pointDensity[rgridNumber] = pointDensity[0];
		pointVelocity[rgridNumber] = pointVelocity[0];
		pointPressure[rgridNumber] = pointPressure[0];
	} else {
	}
}

void Simulation::successiveApproximatonPressure(double& p, double& u, double& R1, double& R2, double& alpha1, double& alpha2, double p1, double p2, double u1, double u2, double rho1, double rho2){
	double c1 = sqrt(gamma*p1/rho1);
	double c2 = sqrt(gamma*p2/rho2);

	if(p > p1){
		alpha1 = sqrt(rho1*((gamma + 1)*p/2 + (gamma - 1)*p1/2));
	} else {
		if(abs(p/p1 - 1) < epsilon){
			alpha1 = rho1*c1;
		} else {
			alpha1 = ((gamma-1)/(2*gamma))*rho1*c1*(1 - p/p1)/(1 - power(p/p1, (gamma-1)/(2*gamma)));
		}
	}
	if(p > p2){
		if(abs(p/p2 - 1) < epsilon){
			alpha2 = rho2*c2;
		} else {
			alpha2 = ((gamma-1)/(2*gamma))*rho2*c2*(1 - p/p2)/(1 - power(p/p2, (gamma-1)/(2*gamma)));
		}
	} else {
		alpha2 = sqrt(rho2*((gamma + 1)*p/2 + (gamma - 1)*p2/2));
	}
	for(int j = 0; j < 25; ++j){
		if(p > p1){
			alpha1 = sqrt(rho1*((gamma + 1)*p/2 + (gamma - 1)*p1/2));
		} else {
			if(abs(p/p1 - 1) < epsilon){
				alpha1 = rho1*c1;
			} else {
				alpha1 = ((gamma-1)/(2*gamma))*rho1*c1*(1 - p/p1)/(1 - power(p/p1, (gamma-1)/(2*gamma)));
			}
		}
		if(p > p2){
			if(abs(p/p2 - 1) < epsilon){
				alpha2 = rho2*c2;
			} else {
				alpha2 = ((gamma-1)/(2*gamma))*rho2*c2*(1 - p/p2)/(1 - power(p/p2, (gamma-1)/(2*gamma)));
			}
		} else {
			alpha2 = sqrt(rho2*((gamma + 1)*p/2 + (gamma - 1)*p2/2));
		}
	}
	u = (alpha1*u1 + alpha2*u2 + p1 - p2)/(alpha1 + alpha2);

	if(p > p1){
		R1 = rho1*((gamma + 1)*p + (gamma - 1)*p1)/((gamma - 1)*p + (gamma + 1)*p1);
	} else {
		R1 = gamma*p/sqr(c1 + (gamma - 1)*(u1 - u)/2);
	}

	if(p > p2){
		R2 = rho1*((gamma + 1)*p + (gamma - 1)*p2)/((gamma - 1)*p + (gamma + 1)*p2);
	} else {
		R2 = gamma*p/sqr(c2 - (gamma - 1)*(u2 - u)/2);
	}
}

void Simulation::TracPen(double* u, double* flux, double cs, double leftFlux, double rightFlux){

}

double Simulation::momentum(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i < rgridNumber) {
		return middleDensity[i]*middleVelocity[i];
	} else {
		printf("i >= rgridNumber");
	}
	return 0;
}

double Simulation::energy(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i < rgridNumber) {
		return middlePressure[i]/(gamma - 1) + middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2;
	} else {
		printf("i >= rgridNumber");
	}
	return 0;
}

double Simulation::temperatureIn(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i < rgridNumber) {
		return middlePressure[i]*massProton/(kBoltzman*middleDensity[i]);
	} else {
		printf("i >= rgridNumber");
	}
	return 0;
}

double Simulation::soundSpeed(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i < rgridNumber) {
		return sqrt(gamma*middlePressure[i]/middleDensity[i]);
	} else {
		printf("i >= rgridNumber");
	}
	return 0;
}

double Simulation::densityFlux(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i <= rgridNumber) {
		return pointDensity[i]*pointVelocity[i];
	} else {
		printf("i > rgridNumber");
	}
	return 0;
}

double Simulation::momentumFlux(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i <= rgridNumber) {
		return pointPressure[i] + pointDensity[i]*pointVelocity[i]*pointVelocity[i];
	} else {
		printf("i > rgridNumber");
	}
	return 0;
}

double Simulation::energyFlux(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i <= rgridNumber) {
		return pointPressure[i]*pointVelocity[i]*gamma/(gamma - 1) + pointDensity[i]*pointVelocity[i]*pointVelocity[i]*pointVelocity[i]/2;
	} else {
		printf("i > rgridNumber");
	}
	return 0;
}

double Simulation::minmod(double a, double b){
	return 0;
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

void Simulation::updateMaxSoundSpeed(){
	maxSoundSpeed = 0;
	for(int i = 0; i < rgridNumber; ++i){
		double cs = (sqrt(gamma*middlePressure[i]/middleDensity[i]) + abs(middleVelocity[i]));
		if(cs > maxSoundSpeed){
			maxSoundSpeed = cs;
		} 
	}
	deltaT = 0.05*deltaR/maxSoundSpeed;
}
