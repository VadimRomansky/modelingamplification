#include "stdafx.h"
#include <list>
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

Simulation::Simulation(){
	initialEnergy = 10E51;
	time = 0;
	tracPen = true;
}

Simulation::~Simulation(){
	delete[] pgrid;
	delete[] grid;
	delete[] middleGrid;
	delete[] deltaR;
	delete[] middleDeltaR;
	delete[] tempGrid;
	delete[] pointDensity;
	delete[] pointVelocity;
	delete[] pointPressure;
	delete[] middleDensity;
	delete[] middleVelocity;
	delete[] middlePressure;

	delete[] tempU;

	for(int i = 0; i < rgridNumber; ++i){
		delete[] distributionFunction[i];
		delete[] tempDistributionFunction[i];
	}
	delete[] distributionFunction;
	delete[] tempDistributionFunction;
}

void Simulation::initializeProfile(){
	downstreamR = 0;
	pgrid = new double[pgridNumber + 1];
	grid = new double[rgridNumber + 1];
	middleGrid = new double[rgridNumber];
	deltaR = new double[rgridNumber];
	middleDeltaR = new double[rgridNumber];
	tempGrid = new double[rgridNumber];
	pointDensity = new double[rgridNumber + 1];
	pointVelocity = new double[rgridNumber + 1];
	pointPressure = new double[rgridNumber + 1];
	middleDensity = new double[rgridNumber];
	middleVelocity = new double[rgridNumber];
	middlePressure = new double[rgridNumber];
	tempU = new double[rgridNumber];
	distributionFunction = new double*[rgridNumber];
	tempDistributionFunction = new double*[rgridNumber];
	for(int i = 0; i < rgridNumber; i ++){
		distributionFunction[i] = new double[pgridNumber];
		tempDistributionFunction[i] = new double[pgridNumber];
		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = 0;
		}
	}
	double r = downstreamR;
	double pressure0 = density0*kBoltzman*temperature/massProton;
	/*switch(simulationType){
	case 1:
		minP = massProton*U0;
		break;
	case 2:
		minP = 2*sqrt(kBoltzman*massProton*temperature);
		break;
	case 3:
		minP = massProton*U0;
		break;
	case 4:
		minP = sqrt(2*massProton*massProton*(initialEnergy/cube(10*deltaR))/density0);
		break;
	default:
		minP = massProton*U0;
		break;
	}*/
	minP = massProton*speed_of_light/1000;
	maxP = minP*10000;

	double tempDeltaR = (upstreamR - downstreamR)/rgridNumber;
	for(int i = 0; i < rgridNumber + 1; ++i){
		grid[i] = r;
		middleGrid[i] = r + tempDeltaR/2;
		tempGrid[i] = grid[i];
		deltaR[i] = tempDeltaR;
		r += tempDeltaR;
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
		case 4 :
			middleDensity[i] = density0;
			middleVelocity[i] = 0;
			if(i <= 10){
				middlePressure[i] = initialEnergy/cube(10*tempDeltaR);
			} else {
				middlePressure[i] = pressure0;
			}
			break;
		default:
			middleDensity[i] = density0;
			if((i > 3*rgridNumber/10) && (i < 5*rgridNumber/10)){
				middleVelocity[i] = U0;
			} else {
				middleVelocity[i] = 0;
			}
			middlePressure[i] = pressure0;
			break;
		}
	}
	double pstep = exp(log(maxP/minP)/pgridNumber);
	pgrid[0] = minP;
	for(int i = 1; i <= pgridNumber; ++i){
		pgrid[i] = pgrid[i - 1]*pstep;
	}
	pgrid[pgridNumber] = maxP;

	for(int i = 0; i < rgridNumber; ++i){
		if(i == 0){
			middleDeltaR[i] = middleGrid[i];
		} else {
			middleDeltaR[i] = middleGrid[i] - middleGrid[i-1];
		}
		for(int j = 0; j < pgridNumber; ++j){
			double p = (pgrid[j] + pgrid[j + 1])/2;
			//distributionFunction[i][j] = (middleDensity[i]/massProton)*exp(-p*p/(2*massProton*kBoltzman*temperatureIn(i)))/cube(sqrt(2*pi*massProton*kBoltzman*temperatureIn(i)));
			//distributionFunction[i][j] = (middleDensity[i]/massProton)*exp(-p*speed_of_light/(kBoltzman*temperatureIn(i)))/cube(sqrt(2*pi*massProton*kBoltzman*temperatureIn(i)));
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
	printf("initialization\n");
	initializeProfile();
	updateMaxSoundSpeed();
	updateShockWavePoint();
	updateParameters();

	printf("creating files\n");
	FILE* outFile;
	FILE* outIteration;
	fopen_s(&outIteration, "./output/iterations.dat","w");
	fclose(outIteration);
	FILE* outExtraIteration;
	fopen_s(&outExtraIteration, "./output/extra_iterations.dat","w");
	fclose(outExtraIteration);
	FILE* outTempGrid;
	fopen_s(&outTempGrid, "./output/temp_grid.dat","w");
	fclose(outTempGrid);
	FILE* outShockWave;
	fopen_s(&outShockWave, "./output/shock_wave.dat","w");
	fclose(outShockWave);
	fopen_s(&outFile, "./output/tamc_radial_profile.dat","w");
	output(outFile,this);
	fclose(outFile);
	FILE* outDistribution;
	FILE* outFullDistribution;
	FILE* outCoordinateDistribution;
	fopen_s(&outDistribution, "./output/distribution.dat","w");
	fopen_s(&outFullDistribution, "./output/fullDistribution.dat","w");
	fopen_s(&outCoordinateDistribution, "./output/coordinateDistribution.dat","w");
	outputDistribution(outDistribution, outFullDistribution, outCoordinateDistribution, this);
	fclose(outCoordinateDistribution);
	fclose(outFullDistribution);
	fclose(outDistribution);
	updateShockWavePoint();
	fopen_s(&outShockWave, "./output/shock_wave.dat","a");
	double shockWaveR = 0;
	if(shockWavePoint >= 0 && shockWavePoint <= rgridNumber){
		shockWaveR = grid[shockWavePoint];
	}

	fprintf(outShockWave, "%d %lf %d %lf\n", 0, time, shockWavePoint, shockWaveR);
	fclose(outShockWave);
	for(int i = 0; i < iterationNumber; ++i){
		printf("iteration ¹ %d\n", i);
		printf("time = %lf\n", time);
		printf("solving\n");
		//evaluateHydrodynamic();
		if(i < 12){
			evaluateHydrodynamic();
		}
		//evaluateCR();
		time = time + deltaT;
		updateShockWavePoint();
		updateGrid();
		updateMaxSoundSpeed();
		updateShockWavePoint();
		updateParameters();
		if(i % writeParameter == 0){
			printf("outputing\n");
			fopen_s(&outFile, "./output/tamc_radial_profile.dat","a");
			output(outFile, this);
			fclose(outFile);

			fopen_s(&outShockWave, "./output/shock_wave.dat","a");
			double shockWaveR = 0;
			if(shockWavePoint >= 0 && shockWavePoint <= rgridNumber){
				shockWaveR = grid[shockWavePoint];
			}

			fprintf(outShockWave, "%d %lf %d %lf\n", i, time, shockWavePoint, shockWaveR);
			fclose(outShockWave);

			fopen_s(&outDistribution, "./output/distribution.dat","w");
			fopen_s(&outFullDistribution, "./output/fullDistribution.dat","a");
			fopen_s(&outCoordinateDistribution, "./output/coordinateDistribution.dat","w");
			outputDistribution(outDistribution, outFullDistribution, outCoordinateDistribution, this);
			fclose(outCoordinateDistribution);
			fclose(outFullDistribution);
			fclose(outDistribution);

			fopen_s(&outIteration, "./output/iterations.dat","a");
			fprintf(outIteration, "%d %28.20lf %28.20lf %28.20lf %28.20lf\n", i, time, mass, totalMomentum, totalEnergy);
			fclose(outIteration);

			fopen_s(&outExtraIteration, "./output/extra_iterations.dat","a");
			fprintf(outExtraIteration, "%d %28.20lf %28.20lf %28.20lf %28.20lf %28.20lf %28.20lf\n", i, time, mass, totalMomentum, totalEnergy, totalKineticEnergy, totalTermalEnergy);
			fclose(outExtraIteration);

			fopen_s(&outTempGrid, "./output/temp_grid.dat","a");
			outputNewGrid(outTempGrid, this);
			fclose(outTempGrid);
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
		//tempDensity[i] = middleDensity[i]*middleGrid[i]*middleGrid[i];
		tempDensity[i] = middleDensity[i];
		tempMomentum[i] = momentum(i)*middleGrid[i]*middleGrid[i];
		//tempMomentum[i] = momentum(i);
		tempEnergy[i] = energy(i)*middleGrid[i]*middleGrid[i];
		//tempEnergy[i] = energy(i);
		dFlux[i] = densityFlux(i);
		mFlux[i] = momentumConvectiveFlux(i);
		eFlux[i] = energyFlux(i);
	}

	TracPen(tempDensity, dFlux, maxSoundSpeed);
	for(int i = 0; i < rgridNumber - 1; ++i){
		tempDensity[i] -= deltaT*2*middleDensity[i]*middleVelocity[i]/middleGrid[i];
	}
	TracPen(tempMomentum, mFlux, maxSoundSpeed);
	for(int i = 0; i < rgridNumber - 1; ++i){
		tempMomentum[i] -= deltaT*middleGrid[i]*middleGrid[i]*(pointPressure[i+1] - pointPressure[i])/(deltaR[i]);
		//tempMomentum[i] -= deltaT*(pointPressure[i+1] - pointPressure[i])/(deltaR);
		//tempMomentum[i] -= deltaT*2*middleDensity[i]*middleVelocity[i]*middleVelocity[i]/middleGrid[i];
	}
	TracPen(tempEnergy, eFlux, maxSoundSpeed);
	/*for(int i = 0; i < rgridNumber - 1; ++i){
		tempEnergy[i] -= deltaT*2*(middlePressure[i]*middleVelocity[i]*gamma/(gamma - 1) + pointDensity[i]*middleVelocity[i]*middleVelocity[i]*middleVelocity[i]/2)/middleGrid[i];
	}*/

	for(int i = 0; i < rgridNumber; ++i){
		//middleDensity[i] = tempDensity[i]/(middleGrid[i]*middleGrid[i]);
		if(i != 0 || tempDensity[i] < middleDensity[i]){
			middleDensity[i] = tempDensity[i];
		}
		alertNaNOrInfinity(middleDensity[i], "density = NaN");
		double middleMomentum = tempMomentum[i]/(middleGrid[i]*middleGrid[i]);
		//double middleMomentum = tempMomentum[i];
		alertNaNOrInfinity(middleMomentum, "momentum = NaN");
		double middleEnergy = tempEnergy[i];
		alertNaNOrInfinity(middleEnergy, "energy = NaN");
		middleVelocity[i] = middleMomentum/middleDensity[i];
		middlePressure[i] = (middleEnergy/(middleGrid[i]*middleGrid[i]) - middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2)*(gamma - 1);
		//middlePressure[i] = (middleEnergy - middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2)*(gamma - 1);
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
	//if(middleVelocity[0] < 0){
		middleVelocity[0] = 0;
	//}

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
	alertNaNOrInfinity(p1,"pressure = NaN");
	alertNaNOrInfinity(p2,"pressure = NaN");
	if(p1 <= p2){
		double c1 = sqrt(gamma*p1/rho1);
		double c2 = sqrt(gamma*p2/rho2);

		p = firstApproximationPressure(rho1, rho2, u1, u2, p1, p2);
		alertNaNOrInfinity(p,"pressure = NaN");

		for(int i = 1; i < 1000; ++i){
			double tempP = p - (pressureFunction(p, p1, rho1) + pressureFunction(p, p2, rho2) - (u1 - u2))/(pressureFunctionDerivative(p, p1, rho1) + pressureFunctionDerivative(p, p2, rho2));
			alertNaNOrInfinity(tempP,"tempPressure = NaN");
			if(tempP < 0){
				p = p/2;
			} else {
				if(abs(tempP/p - 1) < 0.01*epsilon){
					break;
				}
				p = tempP;
			}
		}

		alertNaNOrInfinity(p,"pressure = NaN");
		if(p >= p1){
			alpha1 = sqrt(rho1*((gamma + 1)*p/2 + (gamma - 1)*p1/2));
		} else {
			if(abs(p/p1 - 1) < epsilon){
				alpha1 = rho1*c1;
			} else {
				alpha1 = ((gamma - 1)/(2*gamma))*rho1*c1*((1 - p/p1)/(1 - power(p/p1, (gamma - 1)/(2*gamma))));
			}
		}
		alertNaNOrInfinity(alpha1,"alpha1 = NaN");

		if(p >= p2){
			alpha2 = sqrt(rho2*((gamma + 1)*p/2 + (gamma - 1)*p2/2));
		} else {
			if(abs(p/p2 - 1) < epsilon){
				alpha2 = rho2*c2;
			} else {
				alpha2 = ((gamma - 1)/(2*gamma))*rho2*c2*((1 - p/p2)/(1 - power(p/p2, (gamma - 1)/(2*gamma))));
			}
		}
		alertNaNOrInfinity(alpha2,"alpha2 = NaN");

		u = (alpha1*u1 + alpha2*u2 + p1 - p2)/(alpha1 + alpha2);
		alertNaNOrInfinity(u,"u = NaN");

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
		if(middleDensity[i]*volume(i) - dt*4*pi*(grid[i+1]*grid[i+1]*densityFlux(i+1) - grid[i]*grid[i]*densityFlux(i))< 0){
		//if(middleDensity[i]*volume(i) - dt*4*pi*(densityFlux(i+1) - densityFlux(i))< 0){
			dt = 0.5*(middleDensity[i]*volume(i)/(4*pi*(grid[i+1]*grid[i+1]*densityFlux(i+1) - grid[i]*grid[i]*densityFlux(i))));
			//dt = 0.5*(middleDensity[i]*volume(i)/(4*pi*(densityFlux(i+1) - densityFlux(i))));
			alertNaNOrInfinity(dt, "dt = NaN");
		}
	}
	deltaT = 0.5*dt;
}

void Simulation::TracPen(double* u, double* flux, double cs){
	cs = cs;

	tempU[0] = u[0] - deltaT*(flux[1] - flux[0])/deltaR[0];
	for(int i = 1; i < rgridNumber - 1; ++i){
		tempU[i] = u[i] - deltaT*((flux[i+1] - flux[i]) - cs*(u[i+1] - 2*u[i] + u[i-1])/2)/deltaR[i];
	}
	tempU[rgridNumber - 1] = u[rgridNumber - 1];

	for(int i = 0; i < rgridNumber; ++i){
		u[i] = tempU[i];
	}
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
		//return grid[i]*grid[i]*pointDensity[i]*pointVelocity[i];
		return pointDensity[i]*pointVelocity[i];
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
		//return pointDensity[i]*pointVelocity[i]*pointVelocity[i];
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
		//return pointPressure[i]*pointVelocity[i]*gamma/(gamma - 1) + pointDensity[i]*pointVelocity[i]*pointVelocity[i]*pointVelocity[i]/2;
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
	double pf = pgrid[injectionMomentum];
	//return 0;
	return 0.1*middleDensity[shockWavePoint]*abs(middleVelocity[shockWavePoint])/(pf*pf*massProton);
	//return 0.1*middleDensity[shockWavePoint]*abs(middleVelocity[shockWavePoint])*pf/massProton;
}

void Simulation::evaluateCR(){
	double* upper = new double[rgridNumber];
	double* middle = new double[rgridNumber];
	double* lower = new double[rgridNumber];

	double* f = new double[rgridNumber];
 	double* x = new double[rgridNumber];

	//for(int j = pgridNumber - 1; j > 0; --j){
	for(int j = 0; j < pgridNumber; ++j){
		// -1 ?
		for(int i = 0; i < rgridNumber; ++i){
			double p = (pgrid[j] + pgrid[j+1])/2;
			double deltaP;
			if(j == 0){
				deltaP = pgrid[1] - pgrid[0];
			} else if(j == pgridNumber -1){
				deltaP = pgrid[rgridNumber] - pgrid[pgridNumber - 1];
			} else {
				deltaP = (pgrid[j+1] - pgrid[j-1])/2;
			}
			double volumeFactor = (cube(grid[i+1]) - cube(grid[i]))/3;
			upper[i] = grid[i+1]*grid[i+1]*diffussionCoef(i+1,j)*deltaT/deltaR[i];
			middle[i] = -((grid[i+1]*grid[i+1]*diffussionCoef(i+1,j) + grid[i]*grid[i]*diffussionCoef(i,j))*deltaT/deltaR[i] + volumeFactor);
			lower[i] = grid[i]*grid[i]*diffussionCoef(i,j)*deltaT/deltaR[i];

			double leftDerivative;
			if(j == 0){
				leftDerivative = 0;
			} else {
				leftDerivative = (distributionFunction[i][j] - distributionFunction[i][j-1])*p/(3*deltaP);
			}
			double rightDerivative;
			if(j == pgridNumber - 1){
				rightDerivative = 0;
			} else {
				rightDerivative = (distributionFunction[i][j+1] - distributionFunction[i][j])*p/(3*deltaP);
			}
			double volumeDerivative = (grid[i+1]*grid[i+1]*pointVelocity[i+1] - grid[i]*grid[i]*pointVelocity[i]);

			f[i] = - volumeFactor*distributionFunction[i][j];

			if(volumeDerivative > 0){
				f[i] += - volumeDerivative*rightDerivative*deltaT;
			} else {
				f[i] += - volumeDerivative*leftDerivative*deltaT;
			}

			if(middleVelocity[i] > 0){
				if(i > 0){
					f[i] += volumeFactor*middleVelocity[i]*(distributionFunction[i][j] - distributionFunction[i-1][j])*deltaT/deltaR[i];
					//f[i] += (grid[i+1]*grid[i+1]*pointVelocity[i+1]*distributionFunction[i][j] - grid[i]*grid[i]*pointVelocity[i]*distributionFunction[i-1][j])*deltaT;
					//f[i] += (grid[i+1]*grid[i+1]*middleVelocity[i]*distributionFunction[i][j] - grid[i]*grid[i]*middleVelocity[i-1]*distributionFunction[i-1][j])*deltaT;
				}
			} else {
				if(i < rgridNumber - 1){
					f[i] += volumeFactor*middleVelocity[i]*(distributionFunction[i+1][j] - distributionFunction[i][j])*deltaT/deltaR[i];
					//f[i] += (grid[i+1]*grid[i+1]*pointVelocity[i+1]*distributionFunction[i+1][j] - grid[i]*grid[i]*pointVelocity[i]*distributionFunction[i][j])*deltaT;
					//f[i] += (grid[i+1]*grid[i+1]*middleVelocity[i+1]*distributionFunction[i+1][j] - grid[i]*grid[i]*middleVelocity[i]*distributionFunction[i][j])*deltaT;
				}
			}
			alertNaNOrInfinity(f[i],"f = NaN");
		}
		//middle[rgridNumber - 1] = -((grid[rgridNumber]*grid[rgridNumber]*diffussionCoef(rgridNumber,j) + grid[rgridNumber - 1]*grid[rgridNumber - 1]*diffussionCoef(rgridNumber - 1,j))*deltaT + (cube(grid[rgridNumber]) - cube(grid[rgridNumber - 1]))/3);
		//f[rgridNumber - 1] = -distributionFunction[rgridNumber - 1][j]*(cube(grid[rgridNumber]) - cube(grid[rgridNumber - 1]))/3;
		solveThreeDiagonal(middle, upper, lower, f, x);
		
		for(int i = 0; i < rgridNumber; ++i){
			tempDistributionFunction[i][j] = x[i];
		}
		if((j == injectionMomentum)  && (shockWavePoint >= 0) && (shockWavePoint < rgridNumber)) {
			tempDistributionFunction[shockWavePoint][j] += injection()*deltaT;
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = tempDistributionFunction[i][j];
		}
	}

	delete[] upper;
	delete[] middle;
	delete[] lower;
	delete[] f;
	delete[] x;
}

void Simulation::solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x){
	double* alpha = new double[rgridNumber];
	double* beta = new double[rgridNumber];

	alpha[1] = -upper[0]/middle[0];
	beta[1] = f[0]/middle[0];
	for(int i = 2; i < rgridNumber; ++i){
		alpha[i] = -upper[i-1]/(lower[i]*alpha[i-1] + middle[i-1]);
		beta[i] = (f[i-1] - lower[i]*beta[i-1])/(lower[i]*alpha[i-1] + middle[i-1]);
	}

	x[rgridNumber - 1] = (f[rgridNumber-1] - lower[rgridNumber-2]*beta[rgridNumber-1])/(lower[rgridNumber-2]*alpha[rgridNumber-1] + middle[rgridNumber-1]);
	alertNaNOrInfinity(x[rgridNumber-1],"x = NaN");
	alertNegative(x[rgridNumber-1],"x < 0");

	for(int i = rgridNumber - 2; i >= 0; --i){
		x[i] = alpha[i+1]*x[i+1] + beta[i+1];
		alertNaNOrInfinity(x[i],"x = NaN");
		alertNegative(x[i],"x < 0");
	}

	delete[] alpha;
	delete[] beta;
}

void Simulation::updateMaxSoundSpeed(){
	maxSoundSpeed = (sqrt(gamma*middlePressure[0]/middleDensity[0]) + abs(middleVelocity[0]));
	double tempdt = min(deltaR[0]/maxSoundSpeed, deltaR[1]/maxSoundSpeed);
	double cs = maxSoundSpeed;
	for(int i = 1; i < rgridNumber - 1; ++i){
		cs = (sqrt(gamma*middlePressure[i]/middleDensity[i]) + abs(middleVelocity[i]));
		if(cs > maxSoundSpeed){
			maxSoundSpeed = cs;
		}
		if(deltaR[i]/cs < tempdt){
			tempdt = deltaR[i]/cs;
		}
		if(deltaR[i-1]/cs < tempdt){
			tempdt = deltaR[i-1]/cs;
		}
		if(deltaR[i+1]/cs < tempdt){
			tempdt = deltaR[i+1]/cs;
		}
	}
	cs = (sqrt(gamma*middlePressure[rgridNumber - 1]/middleDensity[rgridNumber - 1]) + abs(middleVelocity[rgridNumber - 1]));
	if(cs > maxSoundSpeed){
		maxSoundSpeed = cs;
	}
	if(deltaR[rgridNumber - 1]/cs < tempdt){
		tempdt = deltaR[rgridNumber - 1]/cs;
	}
	if(deltaR[rgridNumber - 2]/cs < tempdt){
		tempdt = deltaR[rgridNumber - 2]/cs;
	}
	deltaT = 0.1*tempdt;
}

void Simulation::updateShockWavePoint(){
	double maxGrad = epsilon*abs(density0/(grid[rgridNumber] - grid[0]));
	for(int i = 20; i < rgridNumber - 1; ++i){
		double grad = abs((middleDensity[i] - middleDensity[i + 1])/middleDeltaR[i+1]);
		if(grad > maxGrad){
			maxGrad = grad;
			shockWavePoint = i+1;
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

void Simulation::updateGrid(){
	if ((shockWavePoint < 1) || (shockWavePoint > rgridNumber - 1)) return;
	printf("updating grid\n");
	double shockWaveR = grid[shockWavePoint];
	//int rightPoints = -log((1-gridExpLevel)*(grid[rgridNumber] - shockWaveR)/minDeltaR + 1)/log(gridExpLevel)+1;
	//int rightPoints = min(rightPoints, rgridNumber/2);
	double tempGridLevel = 1 + minDeltaR*(grid[rgridNumber]/(rgridNumber*minDeltaR) - 1)/(grid[rgridNumber] - shockWaveR);
	int rightPoints = log(grid[rgridNumber]/(rgridNumber*minDeltaR))/log(tempGridLevel);
	if(rightPoints > rgridNumber/2){
		rightPoints = rgridNumber/2;
		tempGridLevel = findExpLevel((grid[rgridNumber] - shockWaveR)/minDeltaR, rightPoints, 1.0001, 10);
	}
	if(rightPoints <= 1){
		//rightPoints = 2;
		printf("rightPoints <= 1!!!\n");
	}
	int leftPoints = rgridNumber - 1 - rightPoints;

	tempGrid[0] = 0;
	tempGrid[rgridNumber] = grid[rgridNumber];
	//tempGrid[rgridNumber - 1] = grid[rgridNumber] - a0;
	tempGrid[leftPoints] = shockWaveR;
	double leftDeltaR = (shockWaveR - grid[0])/leftPoints;

	for(int i = 1; i < leftPoints; ++i){
		tempGrid[i] = tempGrid[i-1] + leftDeltaR;
	}

	double logLevel = log(tempGridLevel);
	for(int i = leftPoints + 1; i < rgridNumber - 1; ++i){
		tempGrid[i] = tempGrid[i - 1] + minDeltaR*exp((i - leftPoints - 1)*logLevel);
	}
	tempGrid[rgridNumber - 1]= (tempGrid[rgridNumber - 2] + tempGrid[rgridNumber])/2;

	redistributeValues();
}


void Simulation::redistributeValues(){
	int oldCount = 1;
	int newCount = 1;
	double tempDensity = 0;
	double tempMomentum = 0;
	double tempEnergy = 0;
	double* newDensity = new double[rgridNumber];
	double* newMomentum = new double[rgridNumber];
	double* newEnergy = new double[rgridNumber];
	while((oldCount < rgridNumber) || (newCount < rgridNumber)){
		if(tempGrid[newCount] > grid[oldCount]){
			if(tempGrid[newCount - 1] < grid[oldCount - 1]){
				tempDensity += middleDensity[oldCount - 1]*volume(oldCount-1);
				tempMomentum += momentum(oldCount - 1)*volume(oldCount-1);
				tempEnergy += energy(oldCount - 1)*volume(oldCount-1);
			} else {
				tempDensity += middleDensity[oldCount - 1]*4*pi*(cube(grid[oldCount] ) - cube(tempGrid[newCount - 1]))/3;
				tempMomentum += momentum(oldCount - 1)*4*pi*(cube(grid[oldCount] ) - cube(tempGrid[newCount - 1]))/3;
				tempEnergy += energy(oldCount - 1)*4*pi*(cube(grid[oldCount] ) - cube(tempGrid[newCount - 1]))/3;
			}
			++oldCount;
		} else {
			if(tempGrid[newCount - 1] < grid[oldCount - 1]){
				tempDensity += middleDensity[oldCount - 1]*4*pi*(cube(tempGrid[newCount] ) - cube(grid[oldCount - 1]))/3;
				tempMomentum += momentum(oldCount - 1)*4*pi*(cube(tempGrid[newCount] ) - cube(grid[oldCount - 1]))/3;
				tempEnergy += energy(oldCount - 1)*4*pi*(cube(tempGrid[newCount] ) - cube(grid[oldCount - 1]))/3;
			} else {
				tempDensity += middleDensity[oldCount - 1]*4*pi*(cube(tempGrid[newCount] ) - cube(tempGrid[newCount - 1]))/3;
				tempMomentum += momentum(oldCount - 1)*4*pi*(cube(tempGrid[newCount] ) - cube(tempGrid[newCount - 1]))/3;
				tempEnergy += energy(oldCount - 1)*4*pi*(cube(tempGrid[newCount] ) - cube(tempGrid[newCount - 1]))/3;
			}
			newDensity[newCount - 1] = tempDensity;
			alertNaNOrInfinity(newDensity[newCount - 1], "newDensity = NaN");
			alertNegative(newDensity[newCount - 1], "newDensity < 0");
			newMomentum[newCount - 1] = tempMomentum;
			alertNaNOrInfinity(newMomentum[newCount - 1], "newMomentum = NaN");
			newEnergy[newCount - 1] = tempEnergy;
			alertNaNOrInfinity(newEnergy[newCount - 1], "newEnergy = NaN");
			alertNegative(newEnergy[newCount - 1], "newEnergy < 0");
			tempDensity = 0;
			tempMomentum = 0;
			tempEnergy = 0;
			++newCount;
		}
	}
	if(tempGrid[rgridNumber-1] >= grid[rgridNumber-1]){
		tempDensity += middleDensity[rgridNumber - 1]*4*pi*(cube(tempGrid[rgridNumber]) - cube(tempGrid[rgridNumber - 1]))/3;
		tempMomentum += momentum(rgridNumber - 1)*4*pi*(cube(tempGrid[rgridNumber]) - cube(tempGrid[rgridNumber - 1]))/3;
		tempEnergy += energy(rgridNumber - 1)*4*pi*(cube(tempGrid[rgridNumber]) - cube(tempGrid[rgridNumber - 1]))/3;
	}
	newDensity[rgridNumber - 1] = tempDensity;
	alertNaNOrInfinity(newDensity[rgridNumber - 1], "newDensity = NaN");
	alertNegative(newDensity[rgridNumber - 1], "newDensity < 0");
	newMomentum[rgridNumber - 1] = tempMomentum;
	alertNaNOrInfinity(newMomentum[rgridNumber - 1], "newMomentum = NaN");
	newEnergy[rgridNumber - 1] = tempEnergy;
	alertNaNOrInfinity(newEnergy[rgridNumber - 1], "newEnergy = NaN");
	alertNegative(newEnergy[rgridNumber - 1], "newEnergy < 0");

	for(int i = 0; i < rgridNumber; ++i){
		grid[i + 1] = tempGrid[i + 1];
		middleGrid[i] = (grid[i] + grid[i+1])/2;
		deltaR[i] = grid[i+1] - grid[i];
		if(i == 0){
			middleDeltaR[i] = middleGrid[i];
		} else {
			middleDeltaR[i] = middleGrid[i] - middleGrid[i-1];
		}
		newDensity[i] /= volume(i);
		newMomentum[i] /= volume(i);
		newEnergy[i] /= volume(i);

		middleDensity[i] = newDensity[i];
		if(newDensity[i] <= epsilon*density0){
			middleVelocity[i] = 0;
		} else {
			middleVelocity[i] = newMomentum[i]/newDensity[i];
		}
		middlePressure[i] = (newEnergy[i] - middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2)*(gamma - 1);
		alertNegative(middlePressure[i], "middlePressure < 0");
		alertNaNOrInfinity(middlePressure[i], "middlePressure = NaN");
	}

	delete[] newDensity;
	delete[] newMomentum;
	delete[] newEnergy;
}