#include "stdafx.h"
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

Simulation::Simulation(){
	time = 0;
	tracPen = true;
}

Simulation::~Simulation(){
	delete[] pgrid;
	delete[] grid;
	delete[] middleGrid;
	delete[] pointDensity;
	delete[] pointVelocity;
	delete[] pointPressure;
	delete[] middleDensity;
	delete[] middleVelocity;
	delete[] middlePressure;
	for(int i = 0; i < rgridNumber; ++i){
		delete[] distributionFunction[i];
		delete[] tempDistribution[i];
		delete[] distributionDerivative[i];
	}
	delete[] distributionFunction;
	delete[] tempDistribution;
	delete[] distributionDerivative;
}

void Simulation::initializeProfile(){
	downstreamR = 0;
	pgrid = new double[pgridNumber + 1];
	grid = new double[rgridNumber + 1];
	middleGrid = new double[rgridNumber];
	pointDensity = new double[rgridNumber + 1];
	pointVelocity = new double[rgridNumber + 1];
	pointPressure = new double[rgridNumber + 1];
	middleDensity = new double[rgridNumber];
	middleVelocity = new double[rgridNumber];
	middlePressure = new double[rgridNumber];
	distributionFunction = new double*[rgridNumber];
	tempDistribution = new double*[rgridNumber];
	distributionDerivative = new double*[rgridNumber];
	for(int i = 0; i < rgridNumber; i ++){
		distributionFunction[i] = new double[pgridNumber];
		tempDistribution[i] = new double[pgridNumber];
		distributionDerivative[i] = new double[pgridNumber];
		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = 0;
		}
	}
	deltaR = (upstreamR - downstreamR)/rgridNumber;
	double r = downstreamR;
	double pressure0 = density0*kBoltzman*temperature/massProton;
	switch(simulationType){
	case 1:
		minP = 2*massProton*U0;
		shockWavePoint = rgridNumber/10 - 1;
	case 2:
		minP = 2*sqrt(kBoltzman*massProton*temperature);
		shockWavePoint = -1;
	case 3:
		minP = 2*massProton*U0;
		shockWavePoint = rgridNumber/10 - 1;
	default:
		minP = 2*sqrt(kBoltzman*massProton*temperature);
		shockWavePoint = rgridNumber/100 - 1;
	}
	maxP = minP*100;

	for(int i = 0; i < rgridNumber + 1; ++i){
		grid[i] = r;
		middleGrid[i] = r + deltaR/2;
		r += deltaR;
		switch(simulationType){
		case 1 :
			middleDensity[i] = density0;
			if(i < rgridNumber/10){
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
		case 3 :
			middleDensity[i] = density0;
			if(i <= rgridNumber/10){
				middleVelocity[i] = 10*i*U0/rgridNumber;
			} else if(i > rgridNumber/10 && i < 2*rgridNumber/10){
				middleVelocity[i] = U0;
			} else {
				middleVelocity[i] = 0;
			}
			middlePressure[i] = pressure0;
			break;
		default:
			middleDensity[i] = density0;
			middleVelocity[i] = 0;
			if(i <= rgridNumber/1000){
				middlePressure[i] = 100000000*pressure0;
			} else {
				middlePressure[i] = pressure0;
			}
		}
	}
	double pstep = exp(log(maxP/minP)/pgridNumber);
	pgrid[0] = minP;
	for(int i = 1; i <= pgridNumber; ++i){
		pgrid[i] = pgrid[i - 1]*pstep;
	}
	pgrid[pgridNumber] = maxP;

	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < pgridNumber; ++j){
			double p = (pgrid[j] + pgrid[j + 1])/2;
			//distributionFunction[i][j] = (middleDensity[i]/massProton)*exp(-p*speed_of_light/(kBoltzman*temperatureIn(i)));
			//distributionFunction[i][j] = (middleDensity[i]/massProton)*exp(-p*speed_of_light/(kBoltzman*temperatureIn(0)));
			distributionFunction[i][j] = 0;
		}
		pointDensity[i] = middleDensity[i];
		pointVelocity[i] = middleVelocity[i];
		pointPressure[i] = middlePressure[i];
	}

	pointDensity[rgridNumber] = pointDensity[rgridNumber - 1];
	pointVelocity[rgridNumber] = pointVelocity[rgridNumber - 1];
	pointPressure[rgridNumber] = pointPressure[rgridNumber - 1];
	grid[rgridNumber] = upstreamR;
}

void Simulation::simulate(){
	printf("creating files\n");
	FILE* outFile;
	fopen_s(&outFile, "./output/tamc_radial_profile.dat","w");
	fclose(outFile);
	FILE* outIteration;
	fopen_s(&outIteration, "./output/iterations.dat","w");
	fclose(outIteration);
	FILE* outExtraIteration;
	fopen_s(&outExtraIteration, "./output/extra_iterations.dat","w");
	fclose(outExtraIteration);
	FILE* outDistribution;
	fopen_s(&outDistribution, "./output/distribution.dat","w");
	fclose(outDistribution);
	FILE* outFullDistribution;
	fopen_s(&outFullDistribution, "./output/fullDistribution.dat","w");
	fclose(outFullDistribution);
	FILE* outCoordinateDistribution;
	fopen_s(&outCoordinateDistribution, "./output/coordinateDistribution.dat","w");
	fclose(outCoordinateDistribution);
	printf("initialization\n");
	initializeProfile();
	updateMaxSoundSpeed();
	updateParameters();
	output(outFile,this);
	fclose(outFile);
	for(int i = 0; i < iterationNumber; ++i){
		printf("iteration � %d\n", i);
		printf("time = %lf\n", time);
		printf("solving\n");
		evaluateHydrodynamic();
		//evaluateCR();
		time = time + deltaT;
		//updateValues();
		updateMaxSoundSpeed();
		updateShockWavePoint();
		updateParameters();
		if(i % 1000 == 0){
			printf("outputing\n");
			fopen_s(&outFile, "./output/tamc_radial_profile.dat","a");
			output(outFile, this);
			fclose(outFile);
			//fopen_s(&outDistribution, "./output/distribution.dat","w");
			//fopen_s(&outFullDistribution, "./output/fullDistribution.dat","w");
			//fopen_s(&outCoordinateDistribution, "./output/coordinateDistribution.dat","w");
			//outputDistribution(outDistribution, outFullDistribution, outCoordinateDistribution, this);
			//fclose(outCoordinateDistribution);
			//fclose(outFullDistribution);
			//fclose(outDistribution);
			fopen_s(&outIteration, "./output/iterations.dat","a");
			fprintf(outIteration, "%d %28.20lf %28.20lf %28.20lf %28.20lf\n", i, time, mass, totalMomentum, totalEnergy);
			fclose(outIteration);
			fopen_s(&outExtraIteration, "./output/extra_iterations.dat","a");
			fprintf(outExtraIteration, "%d %28.20lf %28.20lf %28.20lf %28.20lf %28.20lf %28.20lf\n", i, time, mass, totalMomentum, totalEnergy, totalKineticEnergy, totalTermalEnergy);
			fclose(outExtraIteration);
		}
	}
}

void Simulation::evaluateHydrodynamic(){
	solveDiscontinious();
	CheckNegativeDensity();

	double* tempDensity = new double[rgridNumber];
	double* tempMomentum = new double[rgridNumber];
	double* tempEnergy = new double[rgridNumber];
	double* dFlux = new double[rgridNumber];
	double* mFlux = new double[rgridNumber];
	double* eFlux = new double[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		tempDensity[i] = middleDensity[i]*middleGrid[i]*middleGrid[i];
		tempMomentum[i] = momentum(i)*middleGrid[i]*middleGrid[i];
		tempEnergy[i] = energy(i)*middleGrid[i]*middleGrid[i];
		dFlux[i] = densityFlux(i);
		mFlux[i] = momentumConvectiveFlux(i);
		eFlux[i] = energyFlux(i);
	}

		TracPen(tempDensity, dFlux, maxSoundSpeed);
		TracPen(tempMomentum, mFlux, maxSoundSpeed);
		//for(int i = 0; i < rgridNumber; ++i){
		for(int i = 0; i < rgridNumber - 1; ++i){
			tempMomentum[i] -= deltaT*middleGrid[i]*middleGrid[i]*(pointPressure[i+1] - pointPressure[i])/(deltaR);
		}
		TracPen(tempEnergy, eFlux, maxSoundSpeed);

		for(int i = 0; i < rgridNumber; ++i){
			middleDensity[i] = tempDensity[i]/(middleGrid[i]*middleGrid[i]);
			alertNaNOrInfinity(middleDensity[i], "density = NaN");
			double middleMomentum = tempMomentum[i]/(middleGrid[i]*middleGrid[i]);
			alertNaNOrInfinity(middleMomentum, "momentum = NaN");
			double middleEnergy = tempEnergy[i];
			alertNaNOrInfinity(middleEnergy, "energy = NaN");
			middleVelocity[i] = middleMomentum/middleDensity[i];
			middlePressure[i] = (middleEnergy/(middleGrid[i]*middleGrid[i]) - middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2)*(gamma - 1);
			if(middleDensity[i] <= epsilon*density0){
				middleDensity[i] = epsilon*density0;
				middleVelocity[i] = 0;
				//middleVelocity[i] = middleMomentum/middleDensity[i];
				middlePressure[i] = middleDensity[i]*kBoltzman*temperature/massProton;
				//middlePressure[i] = (middleEnergy - middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2)*(gamma - 1);
			}
			if(middlePressure[i] <= 0){
				middlePressure[i] = epsilon*middleDensity[i]*kBoltzman*temperature/massProton;
			}
		}
		middleVelocity[0] = 0;

		delete[] tempDensity;
		delete[] tempMomentum;
		delete[] tempEnergy;
		delete[] dFlux;
		delete[] mFlux;
		delete[] eFlux;
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
		double c1 = sqrt(gamma*p1/rho1);
		double c2 = sqrt(gamma*p2/rho2);
		double u;
		double R1;
		double R2;
		double alpha1;
		double alpha2;

		double Uvacuum = -2*c1/(gamma - 1) - 2*c2/(gamma - 1);

		if( u1 - u2 < Uvacuum){
			pointDensity[i] = 0;
			pointVelocity[i] = 0;
			pointPressure[i] = 0;
			continue;
		}

		
		successiveApproximationPressure(p, u, R1, R2, alpha1, alpha2, p1, p2, u1, u2, rho1, rho2);

		bool isLeftShockWave = (p > p1);
		bool isRightShockWave = (p > p2);

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
			pointDensity[i] = middleDensity[i - 1];
			alertNaNOrInfinity(pointDensity[i], "pointDensity = NaN");
			pointVelocity[i] = middleVelocity[i - 1];
			pointPressure[i] = middlePressure[i - 1];
		} else if(D2 < 0){
			pointDensity[i] = middleDensity[i];
			alertNaNOrInfinity(pointDensity[i], "pointDensity = NaN");
			pointVelocity[i] = middleVelocity[i];
			pointPressure[i] = middlePressure[i];
		} else {
			if(u > 0){
				if(isLeftShockWave){
					pointDensity[i] = R1;
					alertNaNOrInfinity(pointDensity[i], "pointDensity = NaN");
					pointVelocity[i] = u;
					pointPressure[i] = p;
				} else {
					double D3 = u - c1 - (gamma - 1)*(u1 - u)/2;
					if(D3 < 0){
						pointDensity[i] = R1;
						alertNaNOrInfinity(pointDensity[i], "pointDensity = NaN");
						pointVelocity[i] = u;
						pointPressure[i] = p;
					} else {
						pointVelocity[i] = (gamma - 1)*middleVelocity[i-1]/(gamma + 1) + 2*c1/(gamma + 1);
						pointPressure[i] = middlePressure[i-1]*power(abs(pointVelocity[i]/c1), 2*gamma/(gamma - 1));
						pointDensity[i] = gamma*pointPressure[i]/(pointVelocity[i]*pointVelocity[i]);
						alertNaNOrInfinity(pointDensity[i], "pointDensity = NaN");
					}
				}
			} else {
				if(isRightShockWave){
					pointDensity[i] = R2;
					alertNaNOrInfinity(pointDensity[i], "pointDensity = NaN");
					pointVelocity[i] = u;
					pointPressure[i] = p;
				} else {
					double D3 = u + c2 - (gamma - 1)*(u2 - u)/2;
					if(D3 < 0){
						pointVelocity[i] = (gamma - 1)*middleVelocity[i]/(gamma + 1) - 2*c2/(gamma + 1);
						pointPressure[i] = middlePressure[i]*power(abs(pointVelocity[i]/c2), 2*gamma/(gamma - 1));
						pointDensity[i] = gamma*pointPressure[i]/(pointVelocity[i]*pointVelocity[i]);
						alertNaNOrInfinity(pointDensity[i], "pointDensity = NaN");
					} else {
						pointDensity[i] = R2;
						alertNaNOrInfinity(pointDensity[i], "pointDensity = NaN");
						pointVelocity[i] = u;
						pointPressure[i] = p;
					}
				}
			}
		}
	}
}

double Simulation::firstApproximationPressure(double rho1, double rho2, double u1, double u2, double p1, double p2){
	double c1 = sqrt(gamma*p1/rho1);
	double c2 = sqrt(gamma*p2/rho2);

	double p = (p1*rho2*c2 + p2*rho1*c1 + (u1 - u2)*rho1*rho2*c1*c2)/(rho1*c1 + rho2*c2);
	if(p > 0){
		return p;
	}
	return (p1 + p2)/2;
}

void Simulation::successiveApproximationPressure(double& p, double& u, double& R1, double& R2, double& alpha1, double& alpha2, double p1, double p2, double u1, double u2, double rho1, double rho2){
	if(p1 <= p2){
		double c1 = sqrt(gamma*p1/rho1);
		double c2 = sqrt(gamma*p2/rho2);

		p = firstApproximationPressure(rho1, rho2, u1, u2, p1, p2);

		for(int i = 1; i < 1000; ++i){
			double tempP = p - (pressureFunction(p, p1, rho1) + pressureFunction(p, p2, rho2) - (u1 - u2))/(pressureFunctionDerivative(p, p1, rho1) + pressureFunctionDerivative(p, p2, rho2));
			if(tempP < 0){
				p = p/2;
			} else {
				if(abs(tempP/p - 1) < 0.01*epsilon){
					break;
				}
				p = tempP;
			}
		}

		if(p >= p1){
			alpha1 = sqrt(rho1*((gamma + 1)*p/2 + (gamma - 1)*p1/2));
		} else {
			if(abs(p/p1 - 1) < epsilon){
				alpha1 = rho1*c1;
			} else {
				alpha1 = ((gamma - 1)/(2*gamma))*rho1*c1*((1 - p/p1)/(1 - power(p/p1, (gamma - 1)/(2*gamma))));
			}
		}

		if(p >= p2){
			alpha2 = sqrt(rho2*((gamma + 1)*p/2 + (gamma - 1)*p2/2));
		} else {
			if(abs(p/p2 - 1) < epsilon){
				alpha2 = rho2*c2;
			} else {
				alpha2 = ((gamma - 1)/(2*gamma))*rho2*c2*((1 - p/p2)/(1 - power(p/p2, (gamma - 1)/(2*gamma))));
			}
		}

		u = (alpha1*u1 + alpha2*u2 + p1 - p2)/(alpha1 + alpha2);

		if(p >= p1){
			R1 = rho1*(((gamma + 1)*p + (gamma - 1)*p1)/((gamma - 1)*p + (gamma + 1)*p1));
		} else {
			R1 = gamma*p/sqr(c1 + (gamma - 1)*(u1 - u)/2);
		}

		if(p >= p2){
			R2 = rho2*(((gamma + 1)*p + (gamma - 1)*p2)/((gamma - 1)*p + (gamma + 1)*p2));
		} else {
			R2 = gamma*p/sqr(c1 - (gamma - 1)*(u2 - u)/2);
		}
		alertNaNOrInfinity(R1, "R1 = NaN\n");
		alertNaNOrInfinity(R2, "R2 = NaN\n");
		alertNaNOrInfinity(p, "point p = NaN\n");
		alertNaNOrInfinity(u, "point u = NaN\n");
	} else {
		successiveApproximationPressure(p, u, R2, R1, alpha2, alpha1, p2, p1, -u2, -u1, rho2, rho1);
		u = -u;
		alpha1 = alpha1;
		alpha2 = alpha2;
	}
}

double Simulation::pressureFunction(double p, double p1, double rho1){
	double c1 = sqrt(gamma*p1/rho1);
	if(p >= p1){
		return (p - p1)/(rho1*c1*sqrt((gamma + 1)*p/(2*gamma*p1) + (gamma - 1)/(2*gamma)));
	} else {
		return 2*c1*(power(p/p1, (gamma - 1)/(2*gamma)) - 1)/(gamma - 1);
	}
}

double Simulation::pressureFunctionDerivative(double p, double p1, double rho1){
	double c1 = sqrt(gamma*p1/rho1);
	if(p >= p1){
		return ((gamma + 1)*p/p1 + 3*gamma - 1)/(4*gamma*rho1*c1*sqrt(cube((gamma + 1)*p/(2*gamma*p1) + (gamma - 1)/(2*gamma))));
	} else {
		return c1*power(p/p1, (gamma - 1)/(2*gamma))/(gamma*p);
	}
}

double Simulation::pressureFunctionDerivative2(double p, double p1, double rho1){
	double c1 = sqrt(gamma*p1/rho1);
	if(p >= p1){
		return - (gamma + 1)*((gamma + 1)*p/p1 + 7*gamma - 1)/(16*gamma*sqr(rho1)*cube(c1)*power((gamma + 1)*p/(2*gamma*p1) + (gamma - 1)/(2*gamma), 2.5));
	} else {
		return - (gamma + 1)*c1*power(p/p1, (gamma - 1)/(2*gamma))/(2*gamma*gamma*p*p);
	}
}

void Simulation::CheckNegativeDensity(){
	double dt = deltaT;
	for(int i = 0; i < rgridNumber; ++i){
		if(middleDensity[i]*volume(i) - dt*(densityFlux(i+1) - densityFlux(i))< 0){
			dt = 0.5*(middleDensity[i]*volume(i)/(densityFlux(i+1) - densityFlux(i)));
			alertNaNOrInfinity(dt, "dt = NaN");
		}
	}
	deltaT = 0.1*dt;
}

void Simulation::TracPen(double* u, double* flux, double cs){
	cs = cs;
    double* uplus = new double[rgridNumber];
	double* uminus = new double[rgridNumber];

	double* fplus = new double[rgridNumber];
	double* fminus = new double[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		uplus[i] = cs*u[i] + flux[i];
		uminus[i] = cs*u[i] - flux[i];
	}

	fplus[0] = 0;
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

	//u[0] -= deltaT*(0.5*(fplus[1] + fminus[1]) - 0.5*(fplus[0] + fminus[0]))/deltaR;
	u[0] -= deltaT*(0.5*(fplus[1] + fminus[1]))/deltaR;
	for(int i = 1; i <= rgridNumber-2; ++i){
		u[i] -= deltaT*0.5*(fplus[i+1] + fminus[i+1] - (fplus[i] + fminus[i]))/deltaR;
	}
	u[rgridNumber - 1] -= deltaT*(- 0.5*(fplus[rgridNumber - 1] + fminus[rgridNumber - 1]))/deltaR;

	delete[] uplus;
	delete[] uminus;

	delete[] fplus;
	delete[] fminus;
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

double Simulation::kineticEnergy(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i < rgridNumber) {
		return middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2;
	} else {
		printf("i >= rgridNumber");
	}
	return 0;
}

double Simulation::termalEnergy(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i < rgridNumber) {
		return middlePressure[i]/(gamma - 1);
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
		if(i == 0) return 0;
		return grid[i]*grid[i]*pointDensity[i]*pointVelocity[i];
	} else {
		printf("i > rgridNumber");
	}
	return 0;
}

double Simulation::momentumConvectiveFlux(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i <= rgridNumber) {
		if(i == 0) return 0;
		return grid[i]*grid[i]*(pointDensity[i]*pointVelocity[i]*pointVelocity[i]);
	} else {
		printf("i > rgridNumber");
	}
	return 0;
}

double Simulation::energyFlux(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i <= rgridNumber) {
		if(i == 0) return 0;
		return grid[i]*grid[i]*(pointPressure[i]*pointVelocity[i]*gamma/(gamma - 1) + pointDensity[i]*pointVelocity[i]*pointVelocity[i]*pointVelocity[i]/2);
	} else {
		printf("i > rgridNumber");
	}
	return 0;
}

double Simulation::volume(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i <= rgridNumber) {
		return 4*pi*(cube(grid[i + 1]) - cube(grid[i]))/3;
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

double Simulation::diffussionCoef(int i, int j){
	double p = pgrid[j];
	double B = B0;
	return p*p/(massProton*electron_charge*B);
}

double Simulation::injection(){
	double pf = pgrid[pgridNumber/10];
	return 0.000000000000001*middleDensity[shockWavePoint]*abs(middleVelocity[shockWavePoint])/(pf*pf*massProton);
}

void Simulation::evaluateCR(){
	double* upper = new double[rgridNumber - 1];
	double* middle = new double[rgridNumber];
	double* lower = new double[rgridNumber - 1];

	double* f = new double[rgridNumber];
	double* x = new double[rgridNumber];

	for(int j = pgridNumber - 1; j > 0; --j){
		// -1 ?
		for(int i = 0; i < rgridNumber - 1; ++i){
			double p = (pgrid[j] + pgrid[j+1])/2;
			double deltaP = (pgrid[j+1] - pgrid[j-1])/2;
			double volumeFactor = (cube(grid[i+1]) - cube(grid[i]))/3;
			upper[i] = grid[i+1]*grid[i+1]*diffussionCoef(i+1,j)/deltaR;
			middle[i] = -((grid[i+1]*grid[i+1]*diffussionCoef(i+1,j) + grid[i]*grid[i]*diffussionCoef(i,j) + deltaR*volumeFactor/deltaT)/deltaR);
			lower[i] = grid[i+1]*grid[i+1]*diffussionCoef(i+1,j)/deltaR;

			f[i] = - volumeFactor*distributionFunction[i][j]/(deltaT) + volumeFactor*middleVelocity[i]*(distributionFunction[i+1][j] -distributionFunction[i][j])/deltaR
				- (distributionFunction[i][j] - distributionFunction[i][j-1])*(p/(3*deltaP))*(grid[i+1]*grid[i+1]*pointVelocity[i+1] - grid[i]*grid[i]*pointVelocity[i]);
			alertNaNOrInfinity(f[i],"f = NaN");
			/*if(abs(f[i]) > epsilon){
				printf("bbb\n");
			}*/
		}
		middle[rgridNumber - 1] = -((grid[rgridNumber]*grid[rgridNumber]*diffussionCoef(rgridNumber,j) + grid[rgridNumber - 1]*grid[rgridNumber - 1]*diffussionCoef(rgridNumber - 1,j) + deltaR*(cube(grid[rgridNumber]) - cube(grid[rgridNumber - 1]))/(3*deltaT))/deltaR);
		f[rgridNumber - 1] = distributionFunction[rgridNumber - 1][j]*(cube(grid[rgridNumber]) - cube(grid[rgridNumber - 1]))/(3*deltaT);
		solveThreeDiagonal(middle, upper, lower, f, x);
		
		for(int i = 0; i < rgridNumber; ++i){
			distributionFunction[i][j] = x[i];
		}
		if((j == 10)  && (shockWavePoint >= 0) && (shockWavePoint < rgridNumber)) {
			distributionFunction[shockWavePoint][j] += injection()*deltaT;
		}
	}

	delete[] upper;
	delete[] middle;
	delete[] lower;
	delete[] f;
	delete[] x;
}

void Simulation::solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x){
	double alpha = f[0]/middle[0];
	double betta = - lower[0]/middle[0];
	x[0] = 0;
	for(int i = 1; i < rgridNumber - 1; ++i){
		x[i] = 0;
		alpha = (f[i] - lower[i-1]*alpha)/(lower[i-1]*betta + middle[i]);
		alertNaNOrInfinity(alpha, "aplpha = NaN");
		betta = - upper[i]/(lower[i-1]*betta + middle[i]);
		alertNaNOrInfinity(betta, "betta = NaN");
	}
	x[rgridNumber - 1] = (f[rgridNumber - 1] - lower[rgridNumber - 2]*alpha)/(lower[rgridNumber-2]*betta + middle[rgridNumber - 1]);
	alertNaNOrInfinity(x[rgridNumber - 1], "x[rgridNumber - 1] = NaN");
	x[rgridNumber - 2] = (f[rgridNumber - 1] - middle[rgridNumber - 1]*x[rgridNumber - 1])/(lower[rgridNumber - 2]);
	alertNaNOrInfinity(x[rgridNumber - 2], "x[rgridNumber - 2] = NaN");
	for(int i = rgridNumber - 3; i >= 0; --i){
		x[i] = (f[i + 1] - middle[i + 1]*x[i + 1] - upper[i + 1]*x[i+2])/lower[i];
		alertNaNOrInfinity(x[i],"x = NaN");
		/*if(abs(x[i]) > epsilon){
			printf("aaaaaa\n");
		}*/
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
	deltaT = 0.5*deltaR/maxSoundSpeed;
}

void Simulation::updateShockWavePoint(){
	double maxGrad = 0;
	for(int i = 9; i < rgridNumber - 1; ++i){
		double grad = middleVelocity[i] - middleVelocity[i + 1];
		if(grad > maxGrad){
			maxGrad = grad;
			shockWavePoint = i;
		}
	}
}

void Simulation::updateParameters(){
	mass = 0;
	totalMomentum = 0;
	totalEnergy = 0;
	totalKineticEnergy = 0;
	totalTermalEnergy = 0;
	for(int i = 0; i < rgridNumber; ++i){
		mass += middleDensity[i]*volume(i);
		totalMomentum += momentum(i)*volume(i);
		totalEnergy += energy(i)*volume(i);
		totalKineticEnergy += kineticEnergy(i)*volume(i);
		totalTermalEnergy += termalEnergy(i)*volume(i);
	}
}