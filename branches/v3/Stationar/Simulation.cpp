#include <list>
#include <time.h>
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

//конструктор
Simulation::Simulation(){
	initialEnergy = 10E49;
	shockWavePoint = -1;
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
	delete[] cosmicRayConcentration;
	delete[] tempDensity;
	delete[] tempMomentum;
	delete[] tempEnergy;
	delete[] vscattering;

	delete[] kgrid;
	delete[] logKgrid;

	delete[] tempU;
	delete[] magneticEnergy;
	delete[] tempMagneticEnergy;

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
	vscattering = new double[rgridNumber];

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
		vscattering[i] = 0;
		dFluxPlus[i] = new double[3];
		dFluxMinus[i] = new double[3];
		mFluxPlus[i] = new double[3];
		mFluxMinus[i] = new double[3];
		eFluxPlus[i] = new double[3];
		eFluxMinus[i] = new double[3];
	}

	cosmicRayPressure = new double[rgridNumber+1];
	cosmicRayConcentration = new double[rgridNumber+1];
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

	double R0 = upstreamR;
	double a = 10000;
	double b = 10000;
	double h1 = 0.5*rgridNumber/log(1.0+a/4);
	double h2 = 0.5*rgridNumber/log(1.0+b/4);
	for(int i = 0; i < rgridNumber/2; ++i){ 
		grid[i] = (2*R0/a)*(1 - exp(-(1.0*(i+1)-0.5*rgridNumber)/h1));
	}
	for(int i=rgridNumber/2; i < rgridNumber; ++i){
		grid[i] = (2*R0/b)*(exp((1.0*(i+1)-0.5*rgridNumber)/h2)-1.0);
	}
	grid[rgridNumber] = upstreamR/2*(1 + 1.0/(rgridNumber));
	for(int i = 1; i <= rgridNumber; ++i){
		if(grid[i] < grid[i-1]){
			printf("grid[i] < grid[i-1]\n");
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		middleGrid[i] = (grid[i] + grid[i+1])/2;
		tempGrid[i] = grid[i];
		deltaR[i] = grid[i+1] - grid[i];
	}
	tempGrid[rgridNumber] = grid[rgridNumber];

	kgrid = new double[kgridNumber];
	logKgrid = new double[kgridNumber];

	minP = 0.5*massProton*speed_of_light;
	maxP = minP*100000000;

	double kmin = max2((1E-9)*electron_charge*B0/(speed_of_light*maxP),0.01/deltaR[rgridNumber/2]);
	double kmax = (1E6)*electron_charge*B0/(speed_of_light*minP);
	if(kmin > kmax/10000000){
		printf("kmin > kmax/10000000\n");
	}
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
	maxRate = new double[rgridNumber];
	magneticEnergy = new double[rgridNumber];
	tempMagneticEnergy = new double[rgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		maxRate[i] = 0;
		magneticField[i] = new double[kgridNumber];
		tempMagneticField[i] = new double[kgridNumber];
		largeScaleField[i] = new double[kgridNumber];
		growth_rate[i] = new double[kgridNumber];
		magneticEnergy[i] = 0;
		for(int k = 0; k < kgridNumber; ++k){
			magneticField[i][k] = (1E-5)*B0*B0*power(1/kgrid[k], 5/3)*power(kgrid[0],2/3);
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
			magneticEnergy[i] += magneticField[i][k]*kgrid[k]*deltaLogK;
		}
		tempMagneticEnergy[i] = magneticEnergy[i];
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
	
	for(int i = 0; i < rgridNumber; ++i){
		double sigma = 4.0;
		double pressure = density0*U0*U0/sigma;
		int count = rgridNumber/2 - 1;
		int intCount = count/10;
		if(i < count){
			middleDensity[i] = density0/sigma;
			middleVelocity[i] = U0;
			middlePressure[i] = pressure*1E-20;
		} else {
			middleDensity[i] = density0;
			middleVelocity[i] = U0/sigma;
			middlePressure[i] = pressure*0.75;
		}
		shockWavePoint = count;
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

	updateDiffusionCoef();
	//test distribution
	/*double Q = (3E-4)*middleDensity[rgridNumber/2-1]*abs(middleVelocity[rgridNumber/2-1]*middleVelocity[rgridNumber/2-1]/speed_of_light)*pgrid[injectionMomentum]/massProton;
	for(int i = 0; i < rgridNumber; ++i){
		double p0 = pgrid[injectionMomentum];
		for(int j = injectionMomentum; j < pgridNumber; ++j){
			double p = pgrid[j];
			if(grid[i] > 0){
				distributionFunction[i][j] = 3*Q/(p*0.75*U0);
			} else{
				distributionFunction[i][j] = 3*Q/(p*0.75*U0)*exp(middleVelocity[i]*grid[i]/diffusionCoef[i][j]);
			}
			tempDistributionFunction[i][j] = distributionFunction[i][j];
		}
	}*/

	pointDensity[rgridNumber] = pointDensity[rgridNumber - 1];
	pointVelocity[rgridNumber] = pointVelocity[rgridNumber - 1];
	pointEnthalpy[rgridNumber] = (energy(rgridNumber-1) + middlePressure[rgridNumber-1])/middleDensity[rgridNumber-1];
	//grid[rgridNumber] = upstreamR;

}

//главная функция
void Simulation::simulate(){
	printf("initialization\n");
	initializeProfile();

	updateDiffusionCoef();

	printf("creating files\n");
	FILE* outIteration = fopen("./output/iterations.dat","w");
	fclose(outIteration);
	FILE* outExtraIteration = fopen("./output/extra_iterations.dat","w");
	fclose(outExtraIteration);
	FILE* outTempGrid = fopen("./output/temp_grid.dat","w");
	fclose(outTempGrid);
	FILE* outShockWave = fopen("./output/shock_wave.dat","w");
	fclose(outShockWave);
	FILE* outFile = fopen("./output/tamc_radial_profile.dat","w");
	output(outFile,this);
	fclose(outFile);
	FILE* outDistribution;
	FILE* outFullDistribution;
	FILE* outCoordinateDistribution;
	outDistribution  = fopen("./output/distribution.dat","w");
	outFullDistribution  = fopen("./output/fullDistribution.dat","w");
	outCoordinateDistribution  = fopen("./output/coordinateDistribution.dat","w");
	outputDistributionP3(outDistribution, outFullDistribution, outCoordinateDistribution, this);
	fclose(outCoordinateDistribution);
	fclose(outFullDistribution);
	fclose(outDistribution);
	FILE* outField = fopen("./output/field.dat","w");
	fclose(outField);
	FILE* coordinateField = fopen("./output/coordinate_field.dat","w");
	fclose(coordinateField);
	FILE* outFullField = fopen("./output/full_field.dat","w");
	fclose(outFullField);
	FILE* outCoef = fopen("./output/diff_coef.dat","w");
	fclose(outCoef);
	FILE* xFile = fopen("./output/xfile.dat","w");
	fclose(xFile);
	FILE* kFile = fopen("./output/kfile.dat","w");
	fclose(kFile);

	currentIteration = 0;
	updateDiffusionCoef();
	//основной цикл
	while(currentIteration < iterationNumber){
		++currentIteration;
		printf("iteration № %d\n", currentIteration);
		//printf("solving\n");

		evaluateHydrodynamic();
		
		evaluateCR();

		evaluateField();
	


		updateAll();

		if(currentIteration % writeParameter == 0){
			//вывод на некоторых итерациях
			printf("outputing\n");
			outFile = fopen("./output/tamc_radial_profile.dat","a");
			output(outFile, this);
			fclose(outFile);

			outDistribution = fopen("./output/distribution.dat","a");
			outFullDistribution = fopen("./output/fullDistribution.dat","a");
			outCoordinateDistribution = fopen("./output/coordinateDistribution.dat","a");

			outputDistributionP3(outDistribution, outFullDistribution, outCoordinateDistribution, this);

			fclose(outCoordinateDistribution);
			fclose(outFullDistribution);
			fclose(outDistribution);

			/*fopen_s(&outDistributionDerivative, "./output/distributionDerivative.dat","a");
			outputDerivativeForDebug(outDistributionDerivative, this);
			fclose(outDistributionDerivative);*/

			outTempGrid = fopen("./output/temp_grid.dat","a");
			outputNewGrid(outTempGrid, this);
			fclose(outTempGrid);

			outFullField = fopen("./output/full_field.dat","w");
			coordinateField = fopen("./output/coordinate_field.dat","a");
			outField = fopen("./output/field.dat","a");
			outCoef = fopen("./output/diff_coef.dat","a");
			xFile = fopen("./output/xfile.dat","w");
			kFile = fopen("./output/kfile.dat","w");
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
	for(int i = 1; i < shockWavePoint; ++i){
		double momentumF = middleDensity[0]*middleVelocity[0]*middleVelocity[0] + middlePressure[0] - (cosmicRayPressure[i] - cosmicRayPressure[0]) - (magneticEnergy[i] - magneticEnergy[0])/2;
		double energyF = (energy(i-1) + middlePressure[i-1])*middleVelocity[i-1] - middleVelocity[i]*(magneticEnergy[i] - magneticEnergy[i-1]);
		for(int k = 0; k < kgridNumber; ++k){
			energyF -= growth_rate[i][k]*magneticField[i][k]*kgrid[k]*deltaLogK*middleDeltaR[i];
		}

		double a = middleDensity[0]*middleVelocity[0]*(0.5 - _gamma/(_gamma-1));
		double b = _gamma*momentumF/(_gamma - 1);
		double c = - energyF;

		double D = b*b - 4*a*c;
		if(D < 0){
			printf("discriminant < 0 i = %d\n",i);
		}
		double u2 = (-b - sqrt(D))/(2*a);
		double u1 = (-b + sqrt(D))/(2*a);

		tempDensity[i] = middleDensity[0]*middleVelocity[0]/u2;
		if(tempDensity[i] < 0){
			printf("density < 0\n");
			tempDensity[i] = epsilon*density0;
		}
		tempMomentum[i] = tempDensity[i]*u2;
		tempEnergy[i] = tempDensity[i]*u2*u2*(0.5 - 1/(_gamma-1)) + momentumF/(_gamma - 1);
		if(tempEnergy[i] < 0){
			printf("energy < 0\n");
			tempEnergy[i] = tempDensity[i]*u2*u2/2;
		}
	}
	double Rtot = middleDensity[shockWavePoint]/tempDensity[shockWavePoint-1];
	double m = tempMomentum[shockWavePoint-1];
	double prevP = (tempEnergy[shockWavePoint-1] - 0.5*sqr(m)/tempDensity[shockWavePoint-1])/(_gamma-1);
	double p = prevP - m*m*((1/middleDensity[shockWavePoint]) - (1/tempDensity[shockWavePoint-1]));
	for(int i = shockWavePoint; i < rgridNumber; ++i){
		tempMomentum[i] = m;
		tempEnergy[i] = (m*m/middleDensity[i]) + p/(_gamma-1);
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
		return middlePressure[i]/(_gamma - 1) + middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2;
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
		return middlePressure[i]/(_gamma - 1);
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
		return sqrt(_gamma*middlePressure[i]/middleDensity[i]);
	} else {
		printf("i >= rgridNumber");
	}
	return 0;
}

void Simulation::updateFluxes(){
	updateFluxes(dFlux, dFluxPlus, dFluxMinus);
	updateFluxes(mFlux, mFluxPlus, mFluxMinus);
	updateFluxes(eFlux, eFluxPlus, eFluxMinus);
}

void Simulation::updateFluxes(double* flux, double** fluxPlus, double** fluxMinus){
	double beta = 0.5;
	double phi = 1.0/3;

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

		middlePressure[i] = (middleEnergy - middleDensity[i]*middleVelocity[i]*middleVelocity[i]/2)*(_gamma - 1);
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

	for(int i = 0; i <= rgridNumber; ++i){
		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = tempDistributionFunction[i][j];
		}
	}
	evaluateCosmicRayPressure();
	evaluateCRFlux();



	//field
	for(int i = 1; i < rgridNumber; ++i){
		for(int k = 0; k < kgridNumber; ++k){
			magneticField[i][k] = tempMagneticField[i][k];
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		magneticEnergy[i] = 0;
		for(int k = 0; k < kgridNumber; ++k){
			magneticEnergy[i] += magneticField[i][k]*kgrid[k]*deltaLogK;
			largeScaleField[i][k] = sqrt(4*pi*magneticEnergy[i] + B0*B0);
		}
		magneticInductionSum[i] = sqrt(4*pi*magneticEnergy[i] + B0*B0);
	}

	updateDiffusionCoef();
	growthRate();
}