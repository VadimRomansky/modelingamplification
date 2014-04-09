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
	injectedParticles = 0;
}

//деструктор
Simulation::~Simulation(){
	delete[] pgrid;
	delete[] logPgrid;
	delete[] grid;
	delete[] middleGrid;
	delete[] deltaR;
	delete[] middleDeltaR;
	delete[] tempGrid;
	delete[] pointDensity;
	delete[] pointVelocity;
	delete[] pointEnthalpy;
	delete[] pointSoundSpeed;
	delete[] middleDensity;
	delete[] middleVelocity;
	delete[] middlePressure;
	delete[] cosmicRayPressure;

	delete[] distrFunDerivative;
	delete[] distrFunDerivative2;

	delete[] tempU;

	for(int i = 0; i < rgridNumber; ++i){
		delete[] distributionFunction[i];
		delete[] tempDistributionFunction[i];
	}
	delete[] distributionFunction;
	delete[] tempDistributionFunction;
}

//инициализация профиля после считывания данных

void Simulation::initializeProfile(){
	downstreamR = 0;
	//upstreamR = 20000;

	distrFunDerivative = new double[pgridNumber];
	distrFunDerivative2 = new double[pgridNumber];

	pgrid = new double[pgridNumber];
	logPgrid = new double[pgridNumber];
	grid = new double[rgridNumber + 1];
	middleGrid = new double[rgridNumber];
	deltaR = new double[rgridNumber];
	middleDeltaR = new double[rgridNumber];
	tempGrid = new double[rgridNumber];
	pointDensity = new double[rgridNumber + 1];
	pointVelocity = new double[rgridNumber + 1];
	pointEnthalpy = new double[rgridNumber + 1];
	pointSoundSpeed = new double[rgridNumber + 1];
	middleDensity = new double[rgridNumber];
	middleVelocity = new double[rgridNumber];
	middlePressure = new double[rgridNumber];
	cosmicRayPressure = new double[rgridNumber+1];
	tempU = new double[rgridNumber];
	distributionFunction = new double*[rgridNumber+1];
	tempDistributionFunction = new double*[rgridNumber+1];
	for(int i = 0; i <= rgridNumber; i ++){
		cosmicRayPressure[i] = 0;
		distributionFunction[i] = new double[pgridNumber];
		tempDistributionFunction[i] = new double[pgridNumber];
		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = 0;
		}
	}
	double r = downstreamR;
	double pressure0 = density0*kBoltzman*temperature/massProton;

	//minP = massProton*speed_of_light/10;
	minP = 0.01*massProton*speed_of_light;
	maxP = minP*10000000;


	deltaR0 = (upstreamR - downstreamR)/rgridNumber;
	grid[0] = -upstreamR/2 + deltaR0;
	for(int i = 1; i <= rgridNumber; ++i){
		grid[i] = grid[i-1] + deltaR0;
	}

	/*double R0 = upstreamR*2/100000;
	double a= upstreamR/2;
	double b = upstreamR/2;
	double h1=0.5*rgridNumber/log(1.0+a/R0);
	double h2=0.5*rgridNumber/log(1.0+b/R0);
	for(int i=0; i < rgridNumber/2; ++ i){ 
		grid[i] = R0*(1 - exp(-(1.0*(i+1)-0.5*rgridNumber)/h1));
	}
	for(int i=rgridNumber/2; i < rgridNumber; ++i){
		grid[i] = R0*(exp((1.0*(i+1)-0.5*rgridNumber)/h1)-1.0);
	}
	grid[rgridNumber] = upstreamR/2*(1 + 1.0/(rgridNumber));*/

	for(int i = 0; i < rgridNumber; ++i){
		middleGrid[i] = (grid[i] + grid[i+1])/2;
		tempGrid[i] = grid[i];
		deltaR[i] = grid[i+1] - grid[i];
	}
	tempGrid[rgridNumber] = grid[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		switch(simulationType){
		//различные варианты профиля
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
				if(i < count){
					middlePressure[i] = initialEnergy/(count*deltaR0);
				} else {
					middlePressure[i] = pressure0;
				}
				shockWavePoint = count;
				shockWaveMoved = true;
			}
			break;
		case 5 :
			{
				double sigma = 4;
				double pressure = density0*U0*U0/sigma;
				int count = rgridNumber/2 - 1;
				int intCount = count/10;
				if(i < count){
					middleDensity[i] = density0/sigma;
					middleVelocity[i] = U0 + 10000000;
					//middleVelocity[i] = 1;
					middlePressure[i] = pressure*0.0000000000001;
				} else {
					middleDensity[i] = density0;
					middleVelocity[i] = U0/sigma + 10000000;
					//middleVelocity[i] = 0.25;
					middlePressure[i] = pressure*0.75;
				}
				/*middleDensity[count - 3] = density0/(0.9*sigma);
				middleVelocity[count -3] = U0*0.9 + 10000000;
				middleDensity[count - 2] = density0/(0.7*sigma);
				middleVelocity[count -2] = U0*0.7 + 10000000;
				middleDensity[count - 1] = density0/(0.5*sigma);
				middleVelocity[count -1] = U0*0.5 + 10000000;*/
				shockWavePoint = count;
				shockWaveMoved = true;
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
	//logPgrid[pgridNumber - 1] = log(maxP);
	deltaLogP = (log(maxP) - logPgrid[0])/(pgridNumber);
	pgrid[0] = minP;
	for(int i = 1; i < pgridNumber; ++i){
		logPgrid[i] = logPgrid[i-1] + deltaLogP;
		pgrid[i] = exp(logPgrid[i]);
	}
	pgrid[pgridNumber-1] = maxP;

	for(int i = 0; i < rgridNumber; ++i){
		if(i == 0){
			middleDeltaR[i] = middleGrid[i];
		} else {
			middleDeltaR[i] = middleGrid[i] - middleGrid[i-1];
		}
		double p;
		for(int j = 0; j < pgridNumber; ++j){
			p = exp(logPgrid[j]);
			distributionFunction[i][j] = 0;
			//double x = -(sqrt(sqr(massProton*speed_of_light*speed_of_light) + sqr(p*speed_of_light))-massProton*speed_of_light*speed_of_light)/(kBoltzman*temperatureIn(i));
			//distributionFunction[i][j] = exp(x);
			//distributionFunction[i][j] = epsilon;
		}
		pointDensity[i] = middleDensity[i];
		pointVelocity[i] = middleVelocity[i];
		pointEnthalpy[i] = middlePressure[i] + energy(i);
	}

	pointDensity[rgridNumber] = pointDensity[rgridNumber - 1];
	pointVelocity[rgridNumber] = pointVelocity[rgridNumber - 1];
	pointEnthalpy[rgridNumber] = pointEnthalpy[rgridNumber - 1];
	//grid[rgridNumber] = upstreamR;
}

//главная функция
void Simulation::simulate(){
	printf("initialization\n");
	initializeProfile();
	//updateShockWavePoint();
	//shockWavePoint = rgridNumber/100;
	//updateGrid();
	updateMaxSoundSpeed();
	updateParameters();
	updateTimeStep();

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
	FILE* outDistributionDerivative;
	fopen_s(&outDistributionDerivative, "./output/distributionDerivative.dat","w");
	fclose(outDistributionDerivative);
	updateShockWavePoint();
	fopen_s(&outShockWave, "./output/shock_wave.dat","a");
	double shockWaveR = 0;
	if(shockWavePoint >= 0 && shockWavePoint <= rgridNumber){
		shockWaveR = grid[shockWavePoint];
	}

	fprintf(outShockWave, "%d %lf %d %lf\n", 0, time, shockWavePoint, shockWaveR);
	fclose(outShockWave);
	deltaT = min2(5000, deltaT);
	//deltaT = 5000;
	//deltaT = 0.001;

	clock_t currentTime = clock();
	clock_t prevTime = currentTime;
	currentIteration = 0;
	//основной цикл
	while(myTime < maxTime && currentIteration < iterationNumber){
		++currentIteration;
		printf("iteration № %d\n", currentIteration);
		printf("time = %lf\n", myTime);
		printf("solving\n");
		deltaT = min2(5000, deltaT);
		//deltaT = 5000;
		//deltaT = 0.001;
		//prevTime = clock();
		evaluateHydrodynamic();
		
		//currentTime = clock();
		//printf("dT evaluating hydro = %lf\n", (currentTime - prevTime)*1.0/CLOCKS_PER_SEC);

		//prevTime = clock();
		//CheckNegativeDistribution();
		//evaluateCR();
		//currentTime = clock();
		//printf("dT evaluating cosmic ray = %lf\n", (currentTime - prevTime)*1.0/CLOCKS_PER_SEC);

		myTime = myTime + deltaT;

		//prevTime = clock();
		//updateGrid();
		//currentTime = clock();
		//printf("dT updating grid = %lf\n", (currentTime - prevTime)*1.0/CLOCKS_PER_SEC);

		updateMaxSoundSpeed();
		updateShockWavePoint();
		updateParameters();
		updateTimeStep();
		if(currentIteration % writeParameter == 0){
			//вывод на некоторых итерациях
			printf("outputing\n");
			fopen_s(&outFile, "./output/tamc_radial_profile.dat","a");
			output(outFile, this);
			fclose(outFile);

			fopen_s(&outShockWave, "./output/shock_wave.dat","a");
			double shockWaveR = 0;
			if(shockWavePoint >= 0 && shockWavePoint <= rgridNumber){
				shockWaveR = grid[shockWavePoint];
			}

			fprintf(outShockWave, "%d %lf %d %lf\n", currentIteration, myTime, shockWavePoint, shockWaveR);
			fclose(outShockWave);

			fopen_s(&outDistribution, "./output/distribution.dat","a");
			fopen_s(&outFullDistribution, "./output/fullDistribution.dat","a");
			fopen_s(&outCoordinateDistribution, "./output/coordinateDistribution.dat","a");

			outputDistributionP3(outDistribution, outFullDistribution, outCoordinateDistribution, this);

			fclose(outCoordinateDistribution);
			fclose(outFullDistribution);
			fclose(outDistribution);

			/*fopen_s(&outDistributionDerivative, "./output/distributionDerivative.dat","a");
			outputDerivativeForDebug(outDistributionDerivative, this);
			fclose(outDistributionDerivative);*/

			fopen_s(&outIteration, "./output/iterations.dat","a");
			fprintf(outIteration, "%d %g %g %g %g %g %g\n", currentIteration, myTime, mass, totalMomentum, totalEnergy, injectedParticles, totalParticles);
			fclose(outIteration);

			fopen_s(&outExtraIteration, "./output/extra_iterations.dat","a");
			fprintf(outExtraIteration, "%d %g %g %g %g %g %g %g %g\n", currentIteration, myTime, mass, totalMomentum, totalEnergy, totalKineticEnergy, totalTermalEnergy, injectedParticles, totalParticles);
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
		double* tflux = flux(i);
		dFlux[i] = tflux[0];
		mFlux[i] = tflux[1];
		eFlux[i] = tflux[2];
		delete[] tflux;
	}

	TracPen(tempDensity, dFlux, 0);

	TracPen(tempMomentum, mFlux, 0);
	for(int i = 0; i < rgridNumber - 1; ++i){
		//tempMomentum[i] -= deltaT*(pointPressure[i+1] - pointPressure[i])/(deltaR[i]);
		//tempMomentum[i] -= deltaT*(cosmicRayPressure[i+1] - cosmicRayPressure[i])/(deltaR[i]);
	}
	TracPen(tempEnergy, eFlux, 0);

	if(tempDensity[rgridNumber - 1] < middleDensity[rgridNumber - 1]){
		printf("aaa\n");
	}
	for(int i = 0; i < rgridNumber; ++i){
		//if(i != 0 || tempDensity[i] < middleDensity[i]){
			middleDensity[i] = tempDensity[i];
		//}
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

	delete[] tempDensity;
	delete[] tempMomentum;
	delete[] tempEnergy;
	delete[] dFlux;
	delete[] mFlux;
	delete[] eFlux;
}


//расчет разрывов, задача Римана

void Simulation::solveDiscontinious(){
	double rho1, rho2, u1, u2, c1, c2;
	for(int i = 1; i < rgridNumber; ++i){
		rho1 = middleDensity[i-1];
		rho2 = middleDensity[i];
		u1 = middleVelocity[i-1];
		u2 = middleVelocity[i];
		c1 = sqrt(gamma*middlePressure[i-1]/rho1);
		c2 = sqrt(gamma*middlePressure[i]/rho2);
		double h1 = (energy(i-1) + middlePressure[i-1])/rho1;
		double h2 = (energy(i) + middlePressure[i])/rho2;
		pointDensity[i] = sqrt(rho1*rho2);
		pointVelocity[i] = (sqrt(rho1)*middleVelocity[i-1] + sqrt(middleDensity[i])*middleVelocity[i])/(sqrt(rho1) + sqrt(rho2));
		pointEnthalpy[i] = (sqr(rho1)*h1 + sqrt(rho2)*h2)/(sqrt(rho1) + sqrt(rho2));
		pointSoundSpeed[i] = sqrt((sqrt(rho1)*c1*c1 + sqrt(rho2)*c2*c2)/(sqrt(rho1) + sqrt(rho2)) + 0.5*(gamma - 1)*(sqrt(rho1*rho2)/sqr(sqrt(rho1)+ sqrt(rho1)))*sqr(u1 - u2));
	}
	pointEnthalpy[0] = (middlePressure[0] + energy(0))/middleDensity[0];
	pointDensity[0] = middleDensity[0];
	//pointVelocity[0] = 0;
	pointVelocity[0] = middleVelocity[0];
	pointSoundSpeed[0] = sqrt(gamma*middlePressure[0]/middleDensity[0]);
}


//проверка шага по времени, чтобы плотность не уменьшалась слишком сильно и не становилась отрицательной

void Simulation::CheckNegativeDensity(){
	double dt = deltaT;
	for(int i = 0; i < rgridNumber; ++i){
		if(middleDensity[i]*volume(i) - dt*4*pi*(densityFlux(i+1) - densityFlux(i))< 0){
			dt = 0.5*(middleDensity[i]*volume(i)/(4*pi*(densityFlux(i+1) - densityFlux(i))));
			alertNaNOrInfinity(dt, "dt = NaN");
		}
	}
	deltaT = dt;
}


//Трак и Пен, немного модифицированные

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

double Simulation::densityFlux(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i <= rgridNumber) {
		if(i == 0) return middleVelocity[0]*middleDensity[0];
		return pointDensity[i]*pointVelocity[i];
		//return middleDensity[i-1]*middleVelocity[i-1];
	} else {
		printf("i > rgridNumber");
	}
	return 0;
}

void Simulation::updateFlux(double* flux){
	double* tempFlux = new double[rgridNumber + 1];
	for(int i = 0; i < rgridNumber-1; ++i){
		double deltaFluxLeft = 0.5*(flux[i+1] - flux[i]);
		double deltaFluxRight = 0.5*(flux[i+2] - flux[i+1]);
		tempFlux[i] = flux[i] + superbee(deltaFluxLeft, deltaFluxRight);
	}
	for(int i = 0; i < rgridNumber-1; ++i){
		flux[i] = tempFlux[i];
	}
	delete[] tempFlux;
}

bool Simulation::CheckShockWave(double& u, double p1, double p2,double u1, double u2, double rho1, double rho2){
	double deltarho = rho1 - rho2;
	if (deltarho == 0) return false;
	u = (rho1*u1 - rho2*u2)/deltarho;
	if ( u != u || 0*u != 0*u){
		printf("NaN\n");
		return false;
	}
	double leftMomentumFlux = rho1*(u1 - u);
	double rightMomentumFlux = rho2*(u2 - u);
	double leftEnergyFlux =(u1-u)*(gamma*p1/(gamma-1) + rho1*(u1-u)*(u1-u)/2);
	double rightEnergyFlux = (u2-u)*(gamma*p2/(gamma-1) + rho2*(u2-u)*(u2-u)/2);
	return (abs(leftMomentumFlux-rightMomentumFlux) < 0.0001*(abs(leftMomentumFlux) + abs(rightMomentumFlux))
			&& abs(leftEnergyFlux-rightEnergyFlux) < 0.0001*(abs(leftEnergyFlux) + abs(rightEnergyFlux)));
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

double* Simulation::flux(int i){
	if(i < 0){
		printf("i < 0");
		return flux(1);
	} else if(i == 0){
		return flux(1);
	} else if(i < rgridNumber){

		double* result = new double[3];

		double leftDflux = middleVelocity[i-1]*middleDensity[i-1];
		double rightDflux = middleVelocity[i]*middleDensity[i];
		double leftMflux = leftDflux*middleVelocity[i-1] + middlePressure[i-1];
		double rightMflux = rightDflux*middleVelocity[i] + middlePressure[i];
		double leftEflux = leftDflux*middleVelocity[i-1]*middleVelocity[i-1]/2 + middleVelocity[i-1]*middlePressure[i-1]*gamma/(gamma-1);
		double rightEflux = rightDflux*middleVelocity[i]*middleVelocity[i]/2 + middleVelocity[i]*middlePressure[i]*gamma/(gamma-1);

		double* deltaS = new double[3];
		double* lambda = new double[3];

		lambda[0] = min2(middleVelocity[i-1] - sqrt(gamma*middlePressure[i-1]/middleDensity[i-1]), pointVelocity[i] - pointSoundSpeed[i]);
		lambda[1] = pointVelocity[i];
		lambda[2] = max2(middleVelocity[i] + sqrt(gamma*middlePressure[i]/middleDensity[i]), pointVelocity[i] + pointSoundSpeed[i]);

		deltaS[0] = ((middlePressure[i] - middlePressure[i-1]) - pointDensity[i]*pointSoundSpeed[i]*(middleVelocity[i] - middleVelocity[i-1]))/(2*sqr(pointSoundSpeed[i]));
		deltaS[1] = (sqr(pointSoundSpeed[i])*(middleDensity[i] - middleDensity[i-1]) - (middlePressure[i] - middlePressure[i-1]))/(2*sqr(pointSoundSpeed[i]));
		deltaS[2] = ((middlePressure[i] - middlePressure[i-1]) + pointDensity[i]*pointSoundSpeed[i]*(middleVelocity[i] - middleVelocity[i-1]))/(2*sqr(pointSoundSpeed[i]));

		double** vectors = new double*[3];
		for(int i = 0; i < 3; ++i){
			vectors[i] = new double[3];
		}

		vectors[0][0] = 1;
		vectors[0][1] = 2;
		vectors[0][2] = 1;
		vectors[1][0] = pointVelocity[i] - pointSoundSpeed[i];
		vectors[1][1] = 2*pointVelocity[i];
		vectors[1][2] = pointVelocity[i] + pointSoundSpeed[i];
		vectors[2][0] = pointEnthalpy[i] - pointVelocity[i]*pointSoundSpeed[i];
		vectors[2][1] = sqr(pointVelocity[i]);
		vectors[2][2] = pointEnthalpy[i] + pointVelocity[i]*pointSoundSpeed[i];

		//result[0] = (leftDflux + rightDflux)/2;
		//result[1] = (leftMflux + rightMflux)/2;
		//result[2] = (leftEflux + rightEflux)/2;

		//result[0] = pointDensity[i]*pointVelocity[i];
		//result[1] = pointDensity[i]*pointVelocity[i]*pointVelocity[i] + pointDensity[i]*pointSoundSpeed[i]*pointSoundSpeed[i]/gamma;
		//result[2] = pointDensity[i]*pointVelocity[i]*pointEnthalpy[i];

		result[0] = middleDensity[i-1]*middleVelocity[i-1];
		result[1] = middleDensity[i-1]*middleVelocity[i-1]*middleVelocity[i-1] + middlePressure[i-1];
		result[2] = middleVelocity[i-1]*(energy(i-1)+ middlePressure[i-1]);

		for(int i = 0; i < 3; ++i){
			for(int j = 0; j < 3; ++j){
				result[i] -= 0.5*abs(lambda[j])*deltaS[j]*vectors[i][j];
			}
		}

		for(int i = 0; i < 3; ++i){
			delete[] vectors[i];
		}
		delete[] vectors;
		delete[] deltaS;
		delete[] lambda;


		return result;
	} else if(i == rgridNumber){
		return flux(rgridNumber - 1);
	} else {
		printf("i > rgridNumber\n");
		return flux(rgridNumber-1);
	}
}

double Simulation::volume(int i){
	if(i < 0){
		printf("i < 0");
	} else if(i >= 0 && i <= rgridNumber) {
		return deltaR[i];
	} else {
		printf("i > rgridNumber");
	}
	return 0;
}

double Simulation::minmod(double a, double b){
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

double Simulation::superbee(double a, double b){
	if(abs(a) >= abs(b)){
		return minmod(a, 2*b);
	} else {
		return minmod(2*a, b);
	}
}

//пересчет шага по времени и максимальной скорости звука

void Simulation::updateMaxSoundSpeed(){
	maxSoundSpeed = (sqrt(gamma*middlePressure[0]/middleDensity[0]) + abs(middleVelocity[0]));
	//double tempdt = min2(deltaR[0]/maxSoundSpeed, deltaR[1]/maxSoundSpeed);
	double cs = maxSoundSpeed;
	for(int i = 1; i < rgridNumber - 1; ++i){
		cs = (sqrt(gamma*middlePressure[i]/middleDensity[i]) + abs(middleVelocity[i]));
		if(cs > maxSoundSpeed){
			maxSoundSpeed = cs;
		}
		/*if(deltaR[i]/cs < tempdt){
			tempdt = deltaR[i]/cs;
		}
		if(deltaR[i-1]/cs < tempdt){
			tempdt = deltaR[i-1]/cs;
		}
		if(deltaR[i+1]/cs < tempdt){
			tempdt = deltaR[i+1]/cs;
		}*/
	}
	cs = (sqrt(gamma*middlePressure[rgridNumber - 1]/middleDensity[rgridNumber - 1]) + abs(middleVelocity[rgridNumber - 1]));
	if(cs > maxSoundSpeed){
		maxSoundSpeed = cs;
	}
	/*if(deltaR[rgridNumber - 1]/cs < tempdt){
		tempdt = deltaR[rgridNumber - 1]/cs;
	}
	if(deltaR[rgridNumber - 2]/cs < tempdt){
		tempdt = deltaR[rgridNumber - 2]/cs;
	}
	deltaT = 0.1*tempdt;*/
}

void Simulation::updateTimeStep(){
	//by dx
	double tempdt = min2(deltaR[0]/maxSoundSpeed, deltaR[1]/maxSoundSpeed);
	for(int i = 1; i < rgridNumber - 1; ++i){
		if(deltaR[i]/maxSoundSpeed < tempdt){
			tempdt = deltaR[i]/maxSoundSpeed;
		}
		if(deltaR[i-1]/maxSoundSpeed < tempdt){
			tempdt = deltaR[i-1]/maxSoundSpeed;
		}
		if(deltaR[i+1]/maxSoundSpeed < tempdt){
			tempdt = deltaR[i+1]/maxSoundSpeed;
		}
	}
	if(deltaR[rgridNumber - 1]/maxSoundSpeed < tempdt){
		tempdt = deltaR[rgridNumber - 1]/maxSoundSpeed;
	}
	if(deltaR[rgridNumber - 2]/maxSoundSpeed < tempdt){
		tempdt = deltaR[rgridNumber - 2]/maxSoundSpeed;
	}

	//by dp
	for(int i = 1; i < rgridNumber; ++i){
		if(tempdt > 0.5*abs(3*deltaLogP*middleDeltaR[i]/(middleVelocity[i] - middleVelocity[i-1]))){
			tempdt = 0.5*abs(3*deltaLogP*middleDeltaR[i]/(middleVelocity[i] - middleVelocity[i-1]));
		}
	}

	deltaT = 0.5*tempdt;
}

//определение точки ударной волны

void Simulation::updateShockWavePoint(){
	int tempShockWavePoint = -1;
	//double maxGrad = density0;
	double maxGrad = U0/upstreamR;
	for(int i = max2(11, shockWavePoint-1); i < 9*rgridNumber/10 - 1; ++i){
		//double grad = abs((middleDensity[i] - middleDensity[i + 1])/middleDeltaR[i+1]);
		double grad = (middleVelocity[i] - middleVelocity[i + 1])/middleDeltaR[i+1];

		//double grad = (middleDensity[i]);
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
	totalParticles = 0;
	for(int i = 0; i < rgridNumber; ++i){
		mass += middleDensity[i]*volume(i);
		totalMomentum += momentum(i)*volume(i);
		totalEnergy += energy(i)*volume(i);
		totalKineticEnergy += kineticEnergy(i)*volume(i);
		totalTermalEnergy += termalEnergy(i)*volume(i);
		for(int j = 0; j < pgridNumber; ++j){
			double dp;
			if(j == 0){
				dp = (pgrid[j + 1] - pgrid[j]);
			} else if(j == pgridNumber -1){
				dp = (pgrid[j] - pgrid[j - 1]);
			} else {
				dp = (pgrid[j + 1] - pgrid[j - 1])/2;
			}
			double dr = 0;
			if(i == 0){
				dr = deltaR[0];
			} else if(i ==rgridNumber-1){
				dr = deltaR[rgridNumber-1];
			} else {
				dr = middleGrid[i] - middleGrid[i-1];
			}
			totalParticles += 4*pi*distributionFunction[i][j]*dr*dp/pgrid[j];
		}
	}
}
