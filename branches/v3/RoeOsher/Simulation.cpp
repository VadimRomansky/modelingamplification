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
	totalEnergy = 0;
	totalKineticEnergy = 0;
	totalTermalEnergy = 0;
	totalMagneticEnergy = 0;
	totalParticleEnergy = 0;
	totalMomentum = 0;
	totalParticles = 0;
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
	delete[] tempDensity;
	delete[] tempMomentum;
	delete[] tempEnergy;

	delete[] kgrid;
	delete[] logKgrid;

	delete[] tempU;

	for(int i = 0; i < rgridNumber; ++i){
		delete[] dFluxPlus[i];
		delete[] dFluxMinus[i];
		delete[] mFluxPlus[i];
		delete[] mFluxMinus[i];
		delete[] eFluxPlus[i];
		delete[] eFluxMinus[i];
	}

	for(int i = 0; i <= rgridNumber; ++i){
		delete[] distributionFunction[i];
		delete[] tempDistributionFunction[i];
		delete[] largeScaleField[i];
		delete[] magneticField[i];
		delete[] tempMagneticField[i];
		delete[] crflux[i];
		delete[] growth_rate[i];
		delete[] diffusionCoef[i];
	}

	delete[] integratedFlux;
	delete[] magneticInductionSum;
	delete[] tempMagneticField;
	delete[] largeScaleField;
	delete[] crflux;
	delete[] growth_rate;
	delete[] magneticField;
	delete[] distributionFunction;
	delete[] tempDistributionFunction;
	delete[] diffusionCoef;
}

//инициализация профиля после считывания данных

void Simulation::initializeProfile(){
	downstreamR = 0;
	//upstreamR = 20000;

	prevShockWavePoint = -1;
	shockWavePoint = -1;

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

	tempDensity = new double[rgridNumber];
	tempMomentum = new double[rgridNumber];
	tempEnergy = new double[rgridNumber];

	dFlux = new double[rgridNumber];
	dFluxPlus = new double*[rgridNumber];
	dFluxMinus = new double*[rgridNumber];

	mFlux = new double[rgridNumber];
	mFluxPlus = new double*[rgridNumber];
	mFluxMinus = new double*[rgridNumber];

	eFlux = new double[rgridNumber];
	eFluxPlus = new double*[rgridNumber];
	eFluxMinus = new double*[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		dFluxPlus[i] = new double[3];
		dFluxMinus[i] = new double[3];
		mFluxPlus[i] = new double[3];
		mFluxMinus[i] = new double[3];
		eFluxPlus[i] = new double[3];
		eFluxMinus[i] = new double[3];
	}

	cosmicRayPressure = new double[rgridNumber+1];
	tempU = new double[rgridNumber];
	integratedFlux = new double[rgridNumber];

	distributionFunction = new double*[rgridNumber+1];
	tempDistributionFunction = new double*[rgridNumber+1];

	for(int i = 0; i <= rgridNumber; i ++){
		cosmicRayPressure[i] = 0;
		distributionFunction[i] = new double[pgridNumber];
		tempDistributionFunction[i] = new double[pgridNumber];
		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = 0;
			tempDistributionFunction[i][j] = 0;
		}
	}

	diffusionCoef = new double*[rgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		diffusionCoef[i] = new double[pgridNumber];
		integratedFlux[i] = 0;
		for(int j =0; j< pgridNumber; ++j){
			diffusionCoef[i][j] = pgrid[j]*speed_of_light*speed_of_light/(electron_charge*B0);
		}
	}

	kgrid = new double[kgridNumber];
	logKgrid = new double[kgridNumber];

	minP = 0.01*massProton*speed_of_light;
	maxP = minP*10000000;

	double kmin = (1E-9)*electron_charge*B0/(speed_of_light*maxP);
	double kmax = (1E6)*electron_charge*B0/(speed_of_light*minP);
	deltaLogK = (log(kmax) - log(kmin))/(kgridNumber - 1);
	for(int i = 0; i < kgridNumber; ++i){
		logKgrid[i] = log(kmin) + i*deltaLogK;
		kgrid[i] = exp(logKgrid[i]);
	}
	magneticField = new double*[rgridNumber];
	tempMagneticField = new double*[rgridNumber];
	largeScaleField = new double*[rgridNumber];
	growth_rate = new double*[rgridNumber];
	magneticInductionSum = new double[rgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		magneticField[i] = new double[kgridNumber];
		tempMagneticField[i] = new double[kgridNumber];
		largeScaleField[i] = new double[kgridNumber];
		growth_rate[i] = new double[kgridNumber];
		for(int k = 0; k < kgridNumber; ++k){
			magneticField[i][k] = (1E-6)*B0*B0*power(1/kgrid[k], 5/3)*power(kgrid[0],2/3);
			if(i >= rgridNumber/2-1){
				magneticField[i][k] *= 8;
			}
			tempMagneticField[i][k] = magneticField[i][k];
			if(k == 0){
				largeScaleField[i][k] = sqrt(4*pi*magneticField[i][k]*kgrid[k]*deltaLogK + B0*B0);
			} else {
				largeScaleField[i][k] = sqrt(4*pi*magneticField[i][k]*kgrid[k]*deltaLogK + sqr(largeScaleField[i][k-1]));
			}
			growth_rate[i][k] = 0;
			alertNaNOrInfinity(magneticField[i][k], "magnetic field = NaN");
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		double magneticEnergy = 0;
		for(int k = 0; k < kgridNumber; ++k){
			magneticEnergy += magneticField[i][k]*kgrid[k]*deltaLogK;
		}
		magneticInductionSum[i] = sqrt(4*pi*magneticEnergy + B0*B0);
	}

	crflux = new double*[rgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		crflux[i] = new double[pgridNumber];
		for(int j = 0;j < pgridNumber; ++j){
			crflux[i][j] = 0;
		}
	}

	double r = downstreamR;
	double pressure0 = density0*kBoltzman*temperature/massProton;

	//minP = massProton*speed_of_light/10;


	//deltaR0 = (upstreamR - downstreamR)/rgridNumber;
	//grid[0] = -upstreamR/2 + deltaR0;
	//for(int i = 1; i <= rgridNumber; ++i){
		//grid[i] = grid[i-1] + deltaR0;
	//}

	double R0 = upstreamR*2/100000;
	double a = upstreamR/2;
	double b = upstreamR/2;
	double h1 = 0.5*rgridNumber/log(1.0+a/R0);
	double h2 = 0.5*rgridNumber/log(1.0+b/R0);
	for(int i = 0; i < rgridNumber/2; ++i){ 
		grid[i] = R0*(1 - exp(-(1.0*(i+1)-0.5*rgridNumber)/h1));
	}
	for(int i=rgridNumber/2; i < rgridNumber; ++i){
		grid[i] = R0*(exp((1.0*(i+1)-0.5*rgridNumber)/h1)-1.0);
	}
	grid[rgridNumber] = upstreamR/2*(1 + 1.0/(rgridNumber));

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
					middleVelocity[i] = U0;
					//middleVelocity[i] = 1;
					middlePressure[i] = pressure*1E-20;
				} else {
					middleDensity[i] = density0;
					middleVelocity[i] = U0/sigma;
					//middleVelocity[i] = 0.25;
					middlePressure[i] = pressure*0.75;
				}
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
		pointEnthalpy[i] = (energy(i) + middlePressure[i])/middleDensity[i];
		tempDensity[i] = middleDensity[i];
		tempMomentum[i] = momentum(i);
		tempEnergy[i] = energy(i);
	}

	pointDensity[rgridNumber] = pointDensity[rgridNumber - 1];
	pointVelocity[rgridNumber] = pointVelocity[rgridNumber - 1];
	pointEnthalpy[rgridNumber] = (energy(rgridNumber-1) + middlePressure[rgridNumber-1])/middleDensity[rgridNumber-1];
	//grid[rgridNumber] = upstreamR;

	shockWaveT = 0;
	shockWaveSpeed = 0;
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
	updateDiffusionCoef();
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
	outputDistributionP3(outDistribution, outFullDistribution, outCoordinateDistribution, this);
	fclose(outCoordinateDistribution);
	fclose(outFullDistribution);
	fclose(outDistribution);
	FILE* outField;
	fopen_s(&outField, "./output/field.dat","w");
	fclose(outField);
	FILE* coordinateField;
	fopen_s(&coordinateField, "./output/coordinate_field.dat","w");
	fclose(coordinateField);
	FILE* outFullField;
	fopen_s(&outFullField, "./output/full_field.dat","w");
	fclose(outFullField);
	FILE* outCoef;
	fopen_s(&outCoef, "./output/diff_coef.dat","w");
	fclose(outCoef);
	FILE* xFile;
	fopen_s(&xFile, "./output/xfile.dat","w");
	fclose(xFile);
	FILE* kFile;
	fopen_s(&kFile, "./output/kfile.dat","w");
	fclose(kFile);
	updateShockWavePoint();
	fopen_s(&outShockWave, "./output/shock_wave.dat","w");
	double shockWaveR = 0;
	if(shockWavePoint >= 0 && shockWavePoint <= rgridNumber){
		shockWaveR = grid[shockWavePoint];
	}

	fprintf(outShockWave, "%d %lf %d %lf\n", 0, myTime, shockWavePoint, shockWaveR, 0, 0);
	fclose(outShockWave);
	deltaT = min2(minT, deltaT);
	//deltaT = 5000;
	//deltaT = 0.001;

	clock_t currentTime = clock();
	clock_t prevTime = currentTime;
	currentIteration = 0;
	updateDiffusionCoef();
	//основной цикл
	while(myTime < maxTime && currentIteration < iterationNumber){
		++currentIteration;
		printf("iteration № %d\n", currentIteration);
		printf("time = %lf\n", myTime);
		printf("solving\n");

		evaluateHydrodynamic();
		
		if(currentIteration > startCRevaluation){
			evaluateCR();
		}

		if(currentIteration > startFieldEvaluation){
			evaluateField();
		}		

		myTime = myTime + deltaT;

		updateAll();

		updateMaxSoundSpeed();
		updateShockWavePoint();
		updateParameters();

		updateTimeStep();
		deltaT = min2(minT, deltaT);
		if(currentIteration % writeParameter == 0){
			//вывод на некоторых итерациях
			printf("outputing\n");
			fopen_s(&outFile, "./output/tamc_radial_profile.dat","a");
			output(outFile, this);
			fclose(outFile);

			fopen_s(&outShockWave, "./output/shock_wave.dat","a");
			double shockWaveR = 0;
			double gasSpeed = 0;
			if(shockWavePoint >= 0 && shockWavePoint <= rgridNumber){
				shockWaveR = grid[shockWavePoint];
				gasSpeed = middleVelocity[shockWavePoint - 1];
			}

			fprintf(outShockWave, "%d %lf %d %lf %lf %lf\n", currentIteration, myTime, shockWavePoint, shockWaveR, shockWaveSpeed, gasSpeed);
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
			fprintf(outExtraIteration, "%d %g %g %g %g %g %g %g %g %g %g\n", currentIteration, myTime, mass, totalMomentum, totalEnergy, totalKineticEnergy, totalTermalEnergy, totalParticleEnergy, totalMagneticEnergy, injectedParticles, totalParticles);
			fclose(outExtraIteration);

			fopen_s(&outTempGrid, "./output/temp_grid.dat","a");
			outputNewGrid(outTempGrid, this);
			fclose(outTempGrid);

			fopen_s(&outFullField, "./output/full_field.dat","w");
			fopen_s(&coordinateField, "./output/coordinate_field.dat","a");
			fopen_s(&outField, "./output/field.dat","a");
			fopen_s(&outCoef, "./output/diff_coef.dat","a");
			fopen_s(&xFile, "./output/xfile.dat","w");
			fopen_s(&kFile, "./output/kfile.dat","w");
			outputField(outField, coordinateField, outFullField, outCoef, xFile, kFile, this);
			fclose(outField);
			fclose(coordinateField);
			fclose(outFullField);
			fclose(outCoef);
			fclose(xFile);
			fclose(kFile);
		}
	}
}

//расчет гидродинамики
void Simulation::evaluateHydrodynamic() {
	printf("evaluating hydrodynamic\n");

	solveDiscontinious();
	CheckNegativeDensity();

	for(int i = 0; i < rgridNumber; ++i){
		tempDensity[i] = middleDensity[i];
		tempMomentum[i] = momentum(i);
		tempEnergy[i] = energy(i);
	}
	evaluateFluxes();

	updateFluxes();

	TracPen(tempDensity, dFlux, 0);

	TracPen(tempMomentum, mFlux, 0);
	/*for(int i = 0; i < rgridNumber - 1; ++i){
		
		//tempMomentum[i] -= deltaT*(cosmicRayPressure[i+1] - cosmicRayPressure[i])/(deltaR[i]);
		
		if(i > 0 && i < rgridNumber-1){
			for(int k = 0; k < kgridNumber; ++k){
				tempMomentum[i] -= deltaT*0.5*(magneticField[i][k] - magneticField[i-1][k])*kgrid[k]*deltaLogK/deltaR[i];
			}
		}
		
	}*/
	TracPen(tempEnergy, eFlux, 0);
	for(int i = 1; i < rgridNumber-1; ++i){
		for(int k = 0; k < kgridNumber; ++k){
			//tempEnergy[i] -= deltaT*0.5*middleVelocity[i]*(magneticField[i][k] - magneticField[i-1][k])*kgrid[k]*deltaLogK/deltaR[i];
		}
	}
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
		pointEnthalpy[i] = (sqrt(rho1)*h1 + sqrt(rho2)*h2)/(sqrt(rho1) + sqrt(rho2));
		pointSoundSpeed[i] = sqrt((sqrt(rho1)*c1*c1 + sqrt(rho2)*c2*c2)/(sqrt(rho1) + sqrt(rho2)) + 0.5*(gamma - 1)*(sqrt(rho1*rho2)/sqr(sqrt(rho1)+ sqrt(rho2)))*sqr(u1 - u2));
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
		if(middleDensity[i]*volume(i) - dt*(densityFlux(i+1) - densityFlux(i))< 0){
			dt = 0.5*(middleDensity[i]*volume(i)/(densityFlux(i+1) - densityFlux(i)));
			alertNaNOrInfinity(dt, "dt = NaN");
		}
	}
	deltaT = dt;
}


//Трак и Пен, немного модифицированные

void Simulation::TracPen(double* u, double* flux, double cs){

	tempU[0] = u[0] - deltaT*(flux[1] - flux[0])/deltaR[0];
	for(int i = 1; i < rgridNumber - 1; ++i){
		tempU[i] = u[i] - deltaT*(flux[i+1] - flux[i])/deltaR[i];
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

void Simulation::evaluateFluxes(){
	dFlux[0] = middleDensity[0]*middleVelocity[0];
	mFlux[0] = dFlux[0]*middleVelocity[0] + middlePressure[0];
	eFlux[0] = middleVelocity[0]*(energy(0) + middlePressure[0]);

	double* deltaS = new double[3];
	double* lambdaPlus = new double[3];
	double* lambdaMinus = new double[3];
	double* lambdaMod = new double[3];
	double** vectors = new double*[3];

	for(int k = 0; k < 3; ++k){
		vectors[k] = new double[3];
	}

	for(int i = 1; i < rgridNumber; ++i){
		double leftDflux = middleVelocity[i-1]*middleDensity[i-1];
		double rightDflux = middleVelocity[i]*middleDensity[i];
		double leftMflux = leftDflux*middleVelocity[i-1] + middlePressure[i-1];
		double rightMflux = rightDflux*middleVelocity[i] + middlePressure[i];
		double leftEflux = leftDflux*middleVelocity[i-1]*middleVelocity[i-1]/2 + middleVelocity[i-1]*middlePressure[i-1]*gamma/(gamma-1);
		double rightEflux = rightDflux*middleVelocity[i]*middleVelocity[i]/2 + middleVelocity[i]*middlePressure[i]*gamma/(gamma-1);

		lambdaMod[0] = min2(middleVelocity[i-1] - sqrt(gamma*middlePressure[i-1]/middleDensity[i-1]), pointVelocity[i] - pointSoundSpeed[i]);
		lambdaMod[1] = pointVelocity[i];
		lambdaMod[2] = max2(middleVelocity[i] + sqrt(gamma*middlePressure[i]/middleDensity[i]), pointVelocity[i] + pointSoundSpeed[i]);

		double lambda1 = pointVelocity[i] - pointSoundSpeed[i];
		double lambda2 = pointVelocity[i];
		double lambda3 = pointVelocity[i] + pointSoundSpeed[i];

		lambdaPlus[0] = lambda1*(lambda1 > 0);
		lambdaPlus[1] = lambda2*(lambda2 > 0);
		lambdaPlus[2] = lambda3*(lambda3 > 0);

		lambdaMinus[0] = lambda1*(lambda1 < 0);
		lambdaMinus[1] = lambda2*(lambda2 < 0);
		lambdaMinus[2] = lambda3*(lambda3 < 0);

		deltaS[0] = ((middlePressure[i] - middlePressure[i-1]) - pointDensity[i]*pointSoundSpeed[i]*(middleVelocity[i] - middleVelocity[i-1]))/(2*sqr(pointSoundSpeed[i]));
		deltaS[1] = (sqr(pointSoundSpeed[i])*(middleDensity[i] - middleDensity[i-1]) - (middlePressure[i] - middlePressure[i-1]))/(2*sqr(pointSoundSpeed[i]));
		deltaS[2] = ((middlePressure[i] - middlePressure[i-1]) + pointDensity[i]*pointSoundSpeed[i]*(middleVelocity[i] - middleVelocity[i-1]))/(2*sqr(pointSoundSpeed[i]));

		vectors[0][0] = 1;
		vectors[0][1] = 2;
		vectors[0][2] = 1;
		vectors[1][0] = pointVelocity[i] - pointSoundSpeed[i];
		vectors[1][1] = 2*pointVelocity[i];
		vectors[1][2] = pointVelocity[i] + pointSoundSpeed[i];
		vectors[2][0] = pointEnthalpy[i] - pointVelocity[i]*pointSoundSpeed[i];
		vectors[2][1] = sqr(pointVelocity[i]);
		vectors[2][2] = pointEnthalpy[i] + pointVelocity[i]*pointSoundSpeed[i];

		dFlux[i] = (leftDflux + rightDflux)/2;
		mFlux[i] = (leftMflux + rightMflux)/2;
		eFlux[i] = (leftEflux + rightEflux)/2;

		//for density
		for(int j = 0; j < 3; ++j){
			dFlux[i] -= 0.5*abs(lambdaMod[j])*deltaS[j]*vectors[0][j];
			dFluxPlus[i][j] = lambdaPlus[j]*deltaS[j]*vectors[0][j];
			dFluxMinus[i][j] = lambdaMinus[j]*deltaS[j]*vectors[0][j];
		}
		//for momentum
		for(int j = 0; j < 3; ++j){
			mFlux[i] -= 0.5*abs(lambdaMod[j])*deltaS[j]*vectors[1][j];
			mFluxPlus[i][j] = lambdaPlus[j]*deltaS[j]*vectors[1][j];
			mFluxMinus[i][j] = lambdaMinus[j]*deltaS[j]*vectors[1][j];
		}
		//for energy
		for(int j = 0; j < 3; ++j){
			eFlux[i] -= 0.5*abs(lambdaMod[j])*deltaS[j]*vectors[2][j];
			eFluxPlus[i][j] = lambdaPlus[j]*deltaS[j]*vectors[2][j];
			eFluxMinus[i][j] = lambdaMinus[j]*deltaS[j]*vectors[2][j];
		}
	}

	for(int k = 0; k < 3; ++k){
		delete[] vectors[k];
	}
	delete[] vectors;
	delete[] deltaS;
	delete[] lambdaPlus;
	delete[] lambdaMinus;
	delete[] lambdaMod;
}

void Simulation::updateFluxes(){
	updateFluxes(dFlux, dFluxPlus, dFluxMinus);
	updateFluxes(mFlux, mFluxPlus, mFluxMinus);
	updateFluxes(eFlux, eFluxPlus, eFluxMinus);
}

void Simulation::updateFluxes(double* flux, double** fluxPlus, double** fluxMinus){
	double beta = 0.5;
	double phi = 0.5;

	for(int i = 2; i < rgridNumber-1; ++i){
		for(int j = 0; j < 3; ++j){
			flux[i] += 0.25*(1+phi)*minmod(fluxPlus[i][j], beta*fluxPlus[i-1][j])
					   +0.25*(1-phi)*minmod(beta*fluxPlus[i-1][j], fluxPlus[i][j])
					   -0.25*(1+phi)*minmod(fluxMinus[i][j], beta*fluxMinus[i+1][j])
					   -0.25*(1-phi)*minmod(beta*fluxMinus[i][j], fluxMinus[i+1][j]);

		}
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
	double cs = maxSoundSpeed;
	for(int i = 1; i < rgridNumber - 1; ++i){
		cs = (sqrt(gamma*middlePressure[i]/middleDensity[i]) + abs(middleVelocity[i]));
		if(cs > maxSoundSpeed){
			maxSoundSpeed = cs;
		}
	}
	cs = (sqrt(gamma*middlePressure[rgridNumber - 1]/middleDensity[rgridNumber - 1]) + abs(middleVelocity[rgridNumber - 1]));
	if(cs > maxSoundSpeed){
		maxSoundSpeed = cs;
	}
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
		if(abs(middleVelocity[i] - middleVelocity[i-1])*tempdt > 0.5*abs(3*deltaLogP*middleDeltaR[i])){
			tempdt = 0.5*abs(3*deltaLogP*middleDeltaR[i]/(middleVelocity[i] - middleVelocity[i-1]));
		}
	}


	double maxDiffusion = 0;
	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < pgridNumber; ++j){
			if(diffusionCoef[i][j] > maxDiffusion){
				maxDiffusion  = diffusionCoef[i][j];
			}
		}
	}

	for(int i = 1; i < rgridNumber; ++i){
		double dx = (grid[i+1] - grid[i-1])/2;
		double dxp=grid[i+1]-grid[i];
		double dxm=grid[i]-grid[i-1];
		double a = ((diffusionCoef[i-1][kgridNumber-1]/dxp + diffusionCoef[i][kgridNumber-1]/dxm) - (diffusionCoef[i-1][kgridNumber-1]/dxm) -(diffusionCoef[i][kgridNumber-1]/dxp))/(2*dx);
		if(1 + tempdt*a < 0){
			tempdt = 0.5/abs(a);
		}
	}


	for(int i = 1; i < rgridNumber-1; ++i){
		double dx = (grid[i+1] - grid[i-1])/2;
		double dxp=grid[i+1]-grid[i];
		double dxm=grid[i]-grid[i-1];
		for(int j = 1; j < pgridNumber; ++j){
			double gkp = distributionFunction[i][j];
			double gkm = distributionFunction[i][j-1];
			double der = (1/(2*dx))*(diffusionCoef[i][j]*(distributionFunction[i+1][j] - distributionFunction[i][j])/dxp - diffusionCoef[i-1][j]*(distributionFunction[i][j] - distributionFunction[i-1][j])/dxm) - (1/dx)*(middleVelocity[i]*distributionFunction[i][j] - middleVelocity[i-1]*distributionFunction[i-1][j]) + (1.0/3)*((middleVelocity[i] - middleVelocity[i-1])/dx)*((gkp - gkm)/deltaLogP);
			if((distributionFunction[i][j] > 0) && (distributionFunction[i][j] + der*tempdt < 0)){
				tempdt = 0.5*abs(distributionFunction[i][j]/der);
				if(tempdt < 0.1){
					printf("ooo\n");
				}
			}
		}
	}

	//by field
	for(int i = 1; i < rgridNumber; ++i){
		for(int k = 0; k < kgridNumber; ++k){
			if(tempdt*growth_rate[i][k] > 2){
				tempdt = 0.5*abs(1/growth_rate[i][k]);
			}
		}
	}
	deltaT = 0.5*tempdt;
}

//определение точки ударной волны

void Simulation::updateShockWavePoint(){
	int tempShockWavePoint = -1;
	shockWaveT += deltaT;
	//double maxGrad = density0;
	double maxGrad = U0/upstreamR;
	for(int i = 10; i < 9*rgridNumber/10 - 1; ++i){
		//double grad = abs((middleDensity[i] - middleDensity[i + 1])/middleDeltaR[i+1]);
		double grad = (middleVelocity[i] - middleVelocity[i + 1])/middleDeltaR[i+1];

		//double grad = (middleDensity[i]);
		if(grad > maxGrad){
			maxGrad = grad;
			tempShockWavePoint = i;
		}
	}
	shockWaveMoved = (tempShockWavePoint != shockWavePoint);
	if(shockWaveMoved && (tempShockWavePoint > -1) && (shockWavePoint > -1)){
		prevShockWavePoint = shockWavePoint;
		shockWaveSpeed = (grid[tempShockWavePoint] - grid[prevShockWavePoint])/(shockWaveT);
		shockWaveT = 0;
	}
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
	totalMagneticEnergy = 0;
	totalParticleEnergy = 0;
	for(int i = 0; i < rgridNumber; ++i){
		mass += middleDensity[i]*volume(i);
		totalMomentum += momentum(i)*volume(i);
		totalKineticEnergy += kineticEnergy(i)*volume(i);
		totalTermalEnergy += termalEnergy(i)*volume(i);
		for(int k = 0; k < kgridNumber; ++k){
			totalMagneticEnergy += magneticField[i][k]*kgrid[k]*deltaLogK*volume(i);
		}
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
			totalParticleEnergy += 4*pi*speed_of_light*distributionFunction[i][j]*dr*dp;
		}
	}
	mass -= deltaT*(middleDensity[0]*middleVelocity[0] - middleDensity[rgridNumber-1]*middleVelocity[rgridNumber-1]);
	totalMomentum -= deltaT*(middleDensity[0]*sqr(middleVelocity[0]) + middlePressure[0] - middleDensity[rgridNumber-1]*sqr(middleVelocity[rgridNumber-1]) - middlePressure[rgridNumber-1]);
	for(int k = 0; k < kgridNumber; ++k){
		totalMagneticEnergy -= deltaT*(magneticField[0][k]*middleVelocity[0] - magneticField[rgridNumber-1][k]*middleVelocity[rgridNumber-1])*kgrid[k]*deltaLogK;
	}
	totalKineticEnergy -= deltaT*(middleDensity[0]*cube(middleVelocity[0]) - middleDensity[rgridNumber-1]*cube(middleVelocity[rgridNumber-1]))/2;
	totalTermalEnergy -= deltaT*(middlePressure[0]*middleVelocity[0] - middlePressure[rgridNumber-1]*middleVelocity[rgridNumber-1])*gamma/(gamma-1);
	totalEnergy = totalTermalEnergy + totalKineticEnergy + totalParticleEnergy + totalMagneticEnergy;
}

void Simulation::updateAll(){
	//hydrodinamic
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

	//cosmic rays

	if(currentIteration > startCRevaluation){
		for(int i = 0; i <= rgridNumber; ++i){
			for(int j = 0; j < pgridNumber; ++j){
				distributionFunction[i][j] = tempDistributionFunction[i][j];
			}
		}
		evaluateCosmicRayPressure();
		evaluateCRFlux();
	}


	//field

	if(currentIteration > startFieldEvaluation){
		for(int i = 1; i < rgridNumber; ++i){
			for(int k = 0; k < kgridNumber; ++k){
				magneticField[i][k] = tempMagneticField[i][k];
			}
		}

		for(int i = 0; i < rgridNumber; ++i){
			double magneticEnergy = 0;
			for(int k = 0; k < kgridNumber; ++k){
				magneticEnergy += magneticField[i][k]*kgrid[k]*deltaLogK;
				largeScaleField[i][k] = sqrt(4*pi*magneticEnergy + B0*B0);
			}
			magneticInductionSum[i] = sqrt(4*pi*magneticEnergy + B0*B0);
		}

		updateDiffusionCoef();
		growthRate();
	}
}