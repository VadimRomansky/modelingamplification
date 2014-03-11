#include <list>
#include <time.h>
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

//конструктор
Simulation::Simulation(){
	initialEnergy = 10E49;
	myTime = 0;
	tracPen = true;
	shockWavePoint = -1;
	shockWaveMoved = false;
}

//деструктор
Simulation::~Simulation(){
	delete[] pgrid;
	delete[] logPgrid;
	delete[] grid;
	delete[] gridsquare;
	delete[] volumeFactor;
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
	delete[] cosmicRayPressure;

	delete[] tempU;

	for(int i = 0; i < rgridNumber; ++i){
		delete[] distributionFunction[i];
		delete[] tempDistributionFunction[i];
	}
	delete[] distributionFunction;
	delete[] tempDistributionFunction;
}

//инициализаци€ профил€ после считывани€ данных

void Simulation::initializeProfile(){
	downstreamR = 0;
	pgrid = new double[pgridNumber];
	logPgrid = new double[pgridNumber];
	grid = new double[rgridNumber + 1];
	gridsquare = new double[rgridNumber + 1];
	volumeFactor = new double[rgridNumber];
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
	cosmicRayPressure = new double[rgridNumber];
	tempU = new double[rgridNumber];
	distributionFunction = new double*[rgridNumber];
	tempDistributionFunction = new double*[rgridNumber];
	for(int i = 0; i < rgridNumber; i ++){
		cosmicRayPressure[i] = 0;
		distributionFunction[i] = new double[pgridNumber];
		tempDistributionFunction[i] = new double[pgridNumber];
		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = 0;
		}
	}
	double r = downstreamR;
	double pressure0 = density0*kBoltzman*temperature/massProton;

	minP = massProton*speed_of_light/10;
	maxP = minP*10000000;

	deltaR0 = (upstreamR - downstreamR)/rgridNumber;
	for(int i = 0; i < rgridNumber + 1; ++i){
		grid[i] = r;
		gridsquare[i] = grid[i]*grid[i];
		middleGrid[i] = r + deltaR0/2;
		tempGrid[i] = grid[i];
		deltaR[i] = deltaR0;
		r += deltaR0;
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
			{
				int count = rgridNumber/100;
				if(i <= count){
					middlePressure[i] = initialEnergy/cube(count*deltaR0);
				} else {
					middlePressure[i] = pressure0;
				}
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
	logPgrid[0] = log(minP);
	logPgrid[pgridNumber - 1] = log(maxP);
	double plogstep = (logPgrid[pgridNumber - 1] - logPgrid[0])/(pgridNumber - 1);
	pgrid[0] = minP;
	for(int i = 1; i < pgridNumber; ++i){
		logPgrid[i] = logPgrid[i-1] + plogstep;
		pgrid[i] = exp(logPgrid[i]);
	}
	pgrid[pgridNumber-1] = maxP;

	for(int i = 0; i < rgridNumber; ++i){
		volumeFactor[i] = (cube(grid[i+1]) - cube(grid[i]))/3;
		if(i == 0){
			middleDeltaR[i] = middleGrid[i];
		} else {
			middleDeltaR[i] = middleGrid[i] - middleGrid[i-1];
		}
		for(int j = 0; j < pgridNumber; ++j){
			double p = exp(logPgrid[j]);
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

//главна€ функци€
void Simulation::simulate(){
	printf("initialization\n");
	initializeProfile();
	updateShockWavePoint();
	shockWavePoint = rgridNumber/10;
	updateGrid();
	updateMaxSoundSpeed();
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
	int i = 0;
	deltaT = min2(5000, deltaT);

	clock_t currentTime = clock();
	clock_t prevTime = currentTime;

	while(myTime < maxTime && i < iterationNumber){
		++i;
		printf("iteration є %d\n", i);
		printf("time = %lf\n", myTime);
		printf("solving\n");
		deltaT = min2(5000, deltaT);
		//prevTime = clock();
		evaluateHydrodynamic();
		//currentTime = clock();
		//printf("dT evaluating hydro = %d\n", currentTime - prevTime);

		//prevTime = clock();
		evaluateCR();
		//currentTime = clock();
		//printf("dT evaluating cosmic ray = %d\n", currentTime - prevTime);

		/*if(i < 20000){
			evaluateHydrodynamic();
		} else {
			evaluateCR();
		}*/

		myTime = myTime + deltaT;

		//prevTime = clock();
		updateGrid();
		//currentTime = clock();
		//printf("dT updating grid = %d\n", currentTime - prevTime);

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

			fprintf(outShockWave, "%d %lf %d %lf\n", i, myTime, shockWavePoint, shockWaveR);
			fclose(outShockWave);

			fopen_s(&outDistribution, "./output/distribution.dat","a");
			fopen_s(&outFullDistribution, "./output/fullDistribution.dat","a");
			fopen_s(&outCoordinateDistribution, "./output/coordinateDistribution.dat","w");
			outputDistribution(outDistribution, outFullDistribution, outCoordinateDistribution, this);
			fclose(outCoordinateDistribution);
			fclose(outFullDistribution);
			fclose(outDistribution);

			fopen_s(&outIteration, "./output/iterations.dat","a");
			fprintf(outIteration, "%d %28.20lf %28.20lf %28.20lf %28.20lf\n", i, myTime, mass, totalMomentum, totalEnergy);
			fclose(outIteration);

			fopen_s(&outExtraIteration, "./output/extra_iterations.dat","a");
			fprintf(outExtraIteration, "%d %28.20lf %28.20lf %28.20lf %28.20lf %28.20lf %28.20lf\n", i, myTime, mass, totalMomentum, totalEnergy, totalKineticEnergy, totalTermalEnergy);
			fclose(outExtraIteration);

			fopen_s(&outTempGrid, "./output/temp_grid.dat","a");
			outputNewGrid(outTempGrid, this);
			fclose(outTempGrid);
		}
	}
}

//расчет гидродинамики
void Simulation::evaluateHydrodynamic() {
	solveDiscontinious();
	CheckNegativeDensity();

	double* tempDensity = new double[rgridNumber];
	double* tempMomentum = new double[rgridNumber];
	double* tempEnergy = new double[rgridNumber];
	double* dFlux = new double[rgridNumber];
	double* mFlux = new double[rgridNumber];
	double* eFlux = new double[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		tempDensity[i] = middleDensity[i];
		tempMomentum[i] = momentum(i);
		tempEnergy[i] = energy(i);
		dFlux[i] = densityFlux(i);
		mFlux[i] = momentumConvectiveFlux(i);
		eFlux[i] = energyFlux(i);
	}

	TracPen(tempDensity, dFlux, maxSoundSpeed);
	for(int i = 0; i < rgridNumber - 1; ++i){
		tempDensity[i] -= deltaT*2*middleDensity[i]*middleVelocity[i]/middleGrid[i];
	}
	TracPenRadial(tempMomentum, mFlux, maxSoundSpeed);
	for(int i = 0; i < rgridNumber - 1; ++i){
		/*if(tempMomentum[i] < 0){
			printf("temp momentum < 0\n");
		}*/
		tempMomentum[i] -= deltaT*(pointPressure[i+1] - pointPressure[i])/(deltaR[i]);
		//left or right?
		tempMomentum[i] -= deltaT*(cosmicRayPressure[i+1] - cosmicRayPressure[i])/(deltaR[i]);
		/*if(tempMomentum[i] < 0){
			printf("temp momentum < 0\n");
		}*/
	}
	TracPenRadial(tempEnergy, eFlux, maxSoundSpeed);

	if(tempDensity[rgridNumber - 1] < middleDensity[rgridNumber - 1]){
		printf("aaa\n");
	}
	for(int i = 0; i < rgridNumber; ++i){
		if(i != 0 || tempDensity[i] < middleDensity[i]){
			middleDensity[i] = tempDensity[i];
		}
		alertNaNOrInfinity(middleDensity[i], "density = NaN");
		double middleMomentum = tempMomentum[i];
		alertNaNOrInfinity(middleMomentum, "momentum = NaN");
		double middleEnergy = tempEnergy[i];
		alertNaNOrInfinity(middleEnergy, "energy = NaN");
		middleVelocity[i] = middleMomentum/middleDensity[i];

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


//расчет разрывов, задача –имана

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
	pointPressure[0] = middlePressure[0];
	pointDensity[0] = middleDensity[0];
	pointVelocity[0] = 0;
}


//начальное приближение дл€ давлени€

double Simulation::firstApproximationPressure(double rho1, double rho2, double u1, double u2, double p1, double p2){
	double c1 = sqrt(gamma*p1/rho1);
	double c2 = sqrt(gamma*p2/rho2);

	double p = (p1*rho2*c2 + p2*rho1*c1 + (u1 - u2)*rho1*rho2*c1*c2)/(rho1*c1 + rho2*c2);
	if(p > 0){
		return p;
	}
	return (p1 + p2)/2;
}


//метод поледовательных приближений дл€ нахождени€ давдени€

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
				alpha2= ((gamma - 1)/(2*gamma))*rho2*c2*((1 - p/p2)/(1 - power(p/p2, (gamma - 1)/(2*gamma))));
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


//вспомогательна€ функци€ дл€ вычислний

double Simulation::pressureFunction(double p, double p1, double rho1){
	double c1 = sqrt(gamma*p1/rho1);
	if(p >= p1){
		return (p - p1)/(rho1*c1*sqrt((gamma + 1)*p/(2*gamma*p1) + (gamma - 1)/(2*gamma)));
	} else {
		return 2*c1*(power(p/p1, (gamma - 1)/(2*gamma)) - 1)/(gamma - 1);
	}
}


// и еЄ производные

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


//проверка шага по времени, чтобы плотность не уменьшалась слишком сильно и не становилась отрицательной

void Simulation::CheckNegativeDensity(){
	double dt = deltaT;
	for(int i = 0; i < rgridNumber; ++i){
		if(middleDensity[i]*volume(i) - dt*4*pi*(gridsquare[i+1]*densityFlux(i+1) - gridsquare[i]*densityFlux(i))< 0){
			dt = 0.5*(middleDensity[i]*volume(i)/(4*pi*(gridsquare[i+1]*densityFlux(i+1) - gridsquare[i]*densityFlux(i))));
			alertNaNOrInfinity(dt, "dt = NaN");
		}
	}
	deltaT = 0.5*dt;
}


//“рак и ѕен, немного модифицированные

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

//с учетом сферичности
void Simulation::TracPenRadial(double* u, double* flux, double cs){
	cs = cs;

	tempU[0] = u[0] - deltaT*(flux[1] - flux[0])/(middleGrid[0]*middleGrid[0]*deltaR[0]);
	for(int i = 1; i < rgridNumber - 1; ++i){
		tempU[i] = u[i] - deltaT*((flux[i+1] - flux[i])/(middleGrid[i]*middleGrid[i]) - cs*(u[i+1] - 2*u[i] + u[i-1])/2)/deltaR[i];
		/*if(tempU[i] < 0){
			printf("temp U < 0\n");
		}*/
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
		return gridsquare[i]*(pointDensity[i]*pointVelocity[i]*pointVelocity[i]);
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
		return gridsquare[i]*(pointPressure[i]*pointVelocity[i]*gamma/(gamma - 1) + pointDensity[i]*pointVelocity[i]*pointVelocity[i]*pointVelocity[i]/2);
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

//инжекционный член
double Simulation::injection(){
	double pf = pgrid[injectionMomentum];
	return middleDensity[shockWavePoint]*abs(middleVelocity[shockWavePoint])/(pf*pf*massProton);
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
		volumeDerivative[i] = (gridsquare[i+1]*pointVelocity[i+1] - gridsquare[i]*pointVelocity[i]);
		dtDivdr[i] = deltaT/deltaR[i];
	}

	double deltaLogP = logPgrid[1] - logPgrid[0];
	
	for(int j = 0; j < pgridNumber; ++j){
		// -1 ?
		for(int i = 0; i < rgridNumber; ++i){
			//double p = (pgrid[j] + pgrid[j+1])/2;
			double p = pgrid[j];

			double diffCoefLeft = diffussionCoef(i,j);
			double diffCoefRight = diffussionCoef(i+1,j);

			upper[i] = gridsquare[i+1]*diffCoefRight*dtDivdr[i];
			//middle[i] = -((gridsquare[i+1]*diffCoefRight + gridsquare[i]*diffCoefLeft)*dtDivdr[i] + volumeFactor[i]);
			lower[i] = gridsquare[i]*diffCoefLeft*dtDivdr[i];
			middle[i] = -volumeFactor[i] - upper[i] - lower[i];

			f[i] = - volumeFactor[i]*distributionFunction[i][j];

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
			f[i] += - volumeDerivative[i]*derivative*deltaT;

			if(middleVelocity[i] > 0){
				if(i > 0){
					f[i] += volumeFactor[i]*middleVelocity[i]*(distributionFunction[i][j] - distributionFunction[i-1][j])*dtDivdr[i];
				}
			} else {
				if(i < rgridNumber - 1){
					f[i] += volumeFactor[i]*middleVelocity[i]*(distributionFunction[i+1][j] - distributionFunction[i][j])*dtDivdr[i];
				}
			}
			alertNaNOrInfinity(f[i],"f = NaN");
		}
		//clock_t prevTime = clock();
		solveThreeDiagonal(middle, upper, lower, f, x, alpha, beta);
		//clock_t currentTime = clock();
		//printf("time solving matrix = %d\n", currentTime - prevTime);
		
		for(int i = 0; i < rgridNumber; ++i){
			tempDistributionFunction[i][j] = x[i];
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = tempDistributionFunction[i][j];
		}
	}

	if(shockWavePoint > 0 && shockWavePoint < rgridNumber){
		distributionFunction[shockWavePoint][injectionMomentum] += injection()*deltaT;
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

//решение трЄх диагональной матрицы
void Simulation::solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x, double* alpha, double* beta){

	alpha[1] = -upper[0]/middle[0];
	beta[1] = f[0]/middle[0];
	for(int i = 2; i < rgridNumber; ++i){
		double temp = lower[i]*alpha[i-1] + middle[i-1];
		alpha[i] = -upper[i-1]/temp;
		beta[i] = (f[i-1] - lower[i]*beta[i-1])/temp;
	}

	x[rgridNumber - 1] = (f[rgridNumber-1] - lower[rgridNumber-2]*beta[rgridNumber-1])/(lower[rgridNumber-2]*alpha[rgridNumber-1] + middle[rgridNumber-1]);
	alertNaNOrInfinity(x[rgridNumber-1],"x = NaN");
	alertNegative(x[rgridNumber-1],"x < 0");

	for(int i = rgridNumber - 2; i >= 0; --i){
		x[i] = alpha[i+1]*x[i+1] + beta[i+1];
		alertNaNOrInfinity(x[i],"x = NaN");
		alertNegative(x[i],"x < 0");
	}
}


//вычисление давлени€ космических лучей

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
		pressure *= 4*pi/3;
		cosmicRayPressure[i] = pressure;
	}

	delete[] partPressure;
}

//пересчет шага по времени и максимальной скорости звука

void Simulation::updateMaxSoundSpeed(){
	maxSoundSpeed = (sqrt(gamma*middlePressure[0]/middleDensity[0]) + abs(middleVelocity[0]));
	double tempdt = min2(deltaR[0]/maxSoundSpeed, deltaR[1]/maxSoundSpeed);
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

//определение точки ударной волны

void Simulation::updateShockWavePoint(){
	int tempShockWavePoint = -1;
	double maxGrad = density0;
	for(int i = max2(11, shockWavePoint-1); i < 9*rgridNumber/10 - 1; ++i){
		//double grad = abs((middleDensity[i] - middleDensity[i + 1])/middleDeltaR[i+1]);
		//double grad = (middleVelocity[i] - middleVelocity[i + 1])/middleDeltaR[i+1];

		double grad = (middleDensity[i]);
		if(grad > maxGrad){
			maxGrad = grad;
			tempShockWavePoint = i+1;
		}
	}
	shockWaveMoved = (tempShockWavePoint != shockWavePoint);
	shockWavePoint = tempShockWavePoint;
}


//подсчет полной массы энергии и импульса

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


//изменение сетки

void Simulation::updateGrid(){
	if ((shockWavePoint < 1) || (shockWavePoint > rgridNumber - 1)) return;
	if( !shockWaveMoved) {
		return;
	}
	printf("updating grid\n");
	double shockWaveR = grid[shockWavePoint];
	double rightR = grid[rgridNumber] - shockWaveR;
	double tempRightGridLevel = (0.5*rightR - minDeltaR)/(0.5*rightR - deltaR0);
	if(tempRightGridLevel < 1){
		return;
	}

	int rightPointsExp = log(deltaR0/minDeltaR)/log(tempRightGridLevel);
	if(rightPointsExp > 6*rgridNumber/10){
		rightPointsExp = 6*rgridNumber/10;

	/*int rightPointsExp = log(deltaR0/minDeltaR)/log(tempRightGridLevel);
	if(rightPointsExp > 7*rgridNumber/10){
		rightPointsExp = 7*rgridNumber/10;*/

	}
	if(rightPointsExp <= 1){
		printf("rightPoints <= 1!!!\n");
	}
	if(rightPointsExp < rgridNumber/10){
		return;
	}
	int rightPointsLinear = rightPointsExp/5;

	//double deltaR0Left = 0.5*shockWaveR/(rgridNumber - rightPointsExp - rightPointsLinear);
	//double tempLeftGridLevel = (0.5*shockWaveR - minDeltaR)/(0.5*shockWaveR - deltaR0Left);
	//int leftPointsExp = log(deltaR0Left/minDeltaR)/log(tempLeftGridLevel);
	//leftPointsExp = min2((rgridNumber - rightPointsExp - rightPointsLinear)/10, leftPointsExp);
	// leftPointsLinear = rgridNumber - leftPointsExp - rightPointsExp - rightPointsLinear;
	int leftPointsLinear = rgridNumber - rightPointsExp - rightPointsLinear;

	tempGrid[0] = 0;
	tempGrid[rgridNumber] = grid[rgridNumber];
	//tempGrid[leftPointsLinear + leftPointsExp] = shockWaveR;
	tempGrid[leftPointsLinear] = shockWaveR;
	//shockWavePoint = leftPointsLinear + leftPointsExp;
	shockWavePoint = leftPointsLinear;

	/*double dRLeft = deltaR0Left*exp(-leftPointsExp*log(tempLeftGridLevel));
	if(dRLeft < deltaR0Left){
		printf("dRLeft < minDeltaR");
	}
	double logLevelLeft = log(tempLeftGridLevel);
	for(int i =  leftPointsLinear + leftPointsExp - 1; i > leftPointsLinear; --i){
		tempGrid[i] = tempGrid[i + 1] - dRLeft*exp(-(i - leftPointsLinear - leftPointsExp + 1)*logLevelLeft);
	}*/

	//double leftDeltaR = (tempGrid[leftPointsLinear + 1] - grid[0])/(leftPointsLinear + 1);

	double leftDeltaR = (shockWaveR - grid[0])/(leftPointsLinear + 1);

	for(int i = 1; i < leftPointsLinear; ++i){
		tempGrid[i] = tempGrid[i-1] + leftDeltaR;
	}
	double dR = deltaR0*exp(-rightPointsExp*log(tempRightGridLevel));
	double logLevel = log(tempRightGridLevel);
	/*for(int i = leftPointsLinear + leftPointsExp + 1; i < leftPointsLinear + leftPointsExp + rightPointsExp; ++i){
		tempGrid[i] = tempGrid[i - 1] + dR*exp((i - leftPointsLinear - leftPointsExp - 1)*logLevel);
	}*/
	//dR = (grid[rgridNumber] - tempGrid[leftPointsLinear + rightPointsExp - 1])/(rightPointsLinear + 1);
	for(int i = leftPointsLinear + 1; i < leftPointsLinear + rightPointsExp; ++i){
		tempGrid[i] = tempGrid[i - 1] + dR*exp((i - leftPointsLinear - 1)*logLevel);
	}
	dR = (grid[rgridNumber] - tempGrid[leftPointsLinear + rightPointsExp - 1])/(rightPointsLinear + 1);
	for(int i = leftPointsLinear + rightPointsExp; i < rgridNumber; ++i){
		tempGrid[i] = tempGrid[i - 1] + dR;
		if(tempGrid[i] > grid[rgridNumber]){
			printf("grid > upstreamR\n");
		}
	}
	tempGrid[rgridNumber - 1] = (tempGrid[rgridNumber - 2] + tempGrid[rgridNumber])/2;
	redistributeValues();
}


//перераспределение величин между €чейками новой сетки

void Simulation::redistributeValues(){
	int oldCount = 1;
	double tempDensity = 0;
	double tempMomentum = 0;
	double tempEnergy = 0;
	double* newDensity = new double[rgridNumber];
	double* newMomentum = new double[rgridNumber];
	double* newEnergy = new double[rgridNumber];
	double** newDistributionFunction = new double*[rgridNumber];
	double* tempFunction = new double[pgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		newDistributionFunction[i] = new double[pgridNumber];
		double newVolume = 4*pi*(cube(tempGrid[i+1]) - cube(tempGrid[i]))/3;
		bool oldChanged = false;
		if(tempGrid[i+1] > grid[oldCount]){
			tempDensity = middleDensity[oldCount - 1]*4*pi*(cube(grid[oldCount]) - cube(tempGrid[i]))/(3*newVolume);
			tempMomentum = momentum(oldCount - 1)*4*pi*(cube(grid[oldCount]) - cube(tempGrid[i]))/(3*newVolume);
			tempEnergy = energy(oldCount - 1)*4*pi*(cube(grid[oldCount]) - cube(tempGrid[i]))/(3*newVolume);
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j]*4*pi*(cube(grid[oldCount]) - cube(tempGrid[i]))/(3*newVolume);
			}
			++oldCount;
			oldChanged = true;
		} else {
			tempDensity = middleDensity[oldCount - 1];
			tempMomentum = momentum(oldCount - 1);
			tempEnergy = energy(oldCount - 1);
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j];
			}
		}
		while(grid[oldCount] < tempGrid[i+1]){
			tempDensity += middleDensity[oldCount - 1]*volume(oldCount - 1)/newVolume;
			tempMomentum += momentum(oldCount - 1)*volume(oldCount - 1)/newVolume;
			tempEnergy += energy(oldCount - 1)*volume(oldCount - 1)/newVolume;
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j]*volume(oldCount - 1)/newVolume;
			}
			++oldCount;
			oldChanged = true;
			if(oldCount > rgridNumber + 1){
				printf("oldCount > rgridNUmber + 1\n");
			}
		}
		if(oldChanged){
			tempDensity += middleDensity[oldCount - 1]*4*pi*(cube(tempGrid[i+1]) - cube(grid[oldCount - 1]))/(3*newVolume);
			tempMomentum += momentum(oldCount - 1)*4*pi*(cube(tempGrid[i+1]) - cube(grid[oldCount - 1]))/(3*newVolume);
			tempEnergy += energy(oldCount - 1)*4*pi*(cube(tempGrid[i+1]) - cube(grid[oldCount - 1]))/(3*newVolume);
			for(int j = 0; j < pgridNumber; j++){
				tempFunction[j] = distributionFunction[oldCount-1][j]*4*pi*(cube(tempGrid[i+1]) - cube(grid[oldCount - 1]))/(3*newVolume);
			}
		}
		newDensity[i] = tempDensity;
		alertNaNOrInfinity(newDensity[i], "newDensity = NaN");
		alertNegative(newDensity[i], "newDensity < 0");
		newMomentum[i] = tempMomentum;
		alertNaNOrInfinity(newMomentum[i], "newMomentum = NaN");
		newEnergy[i] = tempEnergy;
		alertNaNOrInfinity(newEnergy[i], "newEnergy = NaN");
		alertNegative(newEnergy[i], "newEnergy < 0");
		for(int j = 0; j < pgridNumber; ++j){
			newDistributionFunction[i][j] = tempFunction[j];
			alertNaNOrInfinity(newDistributionFunction[i][j], "newDistribution = NaN");
			alertNegative(newDistributionFunction[i][j], "newDistribution < 0");
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		grid[i + 1] = tempGrid[i + 1];
		gridsquare[i+1] = grid[i+1]*grid[i+1];
		middleGrid[i] = (grid[i] + grid[i+1])/2;
		volumeFactor[i] = (cube(grid[i+1]) - cube(grid[i]))/3;
		deltaR[i] = grid[i+1] - grid[i];
		if(i == 0){
			middleDeltaR[i] = middleGrid[i];
		} else {
			middleDeltaR[i] = middleGrid[i] - middleGrid[i-1];
		}

		middleDensity[i] = newDensity[i];
		if(newDensity[i] <= epsilon*density0){
			middleVelocity[i] = 0;
		} else {
			middleVelocity[i] = newMomentum[i]/newDensity[i];
		}
		double tempPressure = (newEnergy[i] - middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2)*(gamma - 1);
		
		if(tempPressure < 0){
			middlePressure[i] = 0.01*min2(middlePressure[i+1],middlePressure[i]);
			printf("pressure < 0\n");
		} else {
			middlePressure[i] = tempPressure;
		}
		alertNegative(middlePressure[i], "middlePressure < 0");
		alertNaNOrInfinity(middlePressure[i], "middlePressure = NaN");

		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = newDistributionFunction[i][j];
		}
	}
	if(newDensity[rgridNumber - 1] < middleDensity[rgridNumber - 1]){
		printf("aaa\n");
	}

	for(int i = 0; i < rgridNumber; ++i){
		delete[] newDistributionFunction[i];
	}
	delete[] newDistributionFunction;
	delete[] newDensity;
	delete[] newMomentum;
	delete[] newEnergy;
	delete[] tempFunction;
}
