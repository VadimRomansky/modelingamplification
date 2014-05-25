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
	totalEnergy = 0;
	totalKineticEnergy = 0;
	totalTermalEnergy = 0;
	totalMomentum = 0;
}

//деструктор
Simulation::~Simulation(){
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
	delete[] tempDensity;
	delete[] tempMomentum;
	delete[] tempEnergy;
	delete[] tempU;

	for(int i = 0; i < rgridNumber; ++i){
		delete[] dFluxPlus[i];
		delete[] dFluxMinus[i];
		delete[] mFluxPlus[i];
		delete[] mFluxMinus[i];
		delete[] eFluxPlus[i];
		delete[] eFluxMinus[i];
	}
}

//инициализаци€ профил€ после считывани€ данных

void Simulation::initializeProfile(){
	downstreamR = 0;
	//upstreamR = 20000;

	prevShockWavePoint = -1;
	shockWavePoint = -1;
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

	tempU = new double[rgridNumber];

	/*double R0 = upstreamR;
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
	}*/

	double dR = upstreamR/rgridNumber;
	grid[0] = -upstreamR/2;
	for(int i = 1; i <= rgridNumber; ++i){
		grid[i] = grid[i-1] + dR;
	}

	for(int i = 0; i < rgridNumber; ++i){
		middleGrid[i] = (grid[i] + grid[i+1])/2;
		tempGrid[i] = grid[i];
		deltaR[i] = grid[i+1] - grid[i];
	}
	tempGrid[rgridNumber] = grid[rgridNumber];
	
	int count = rgridNumber/4 - 1;
	for(int i = 0; i < rgridNumber; ++i){
		if(i < count){
			middleDensity[i] = 2;
			middleVelocity[i] = 4;
			middlePressure[i] = 0.0000001;
		} else {
			middleDensity[i] = 0.25;
			middleVelocity[i] = 1;
			middlePressure[i] = 0.75;
		}
		shockWavePoint = count;
		shockWaveMoved = true;
	}

	for(int i = 0; i < rgridNumber; ++i){
		if(i == 0){
			middleDeltaR[i] = middleGrid[i];
		} else {
			middleDeltaR[i] = middleGrid[i] - middleGrid[i-1];
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

	rCD=grid[shockWavePoint];
	rL1=grid[shockWavePoint];
	rL2=grid[shockWavePoint];
	rR1=grid[shockWavePoint];
	rR2=grid[shockWavePoint];
}

//главна€ функци€
void Simulation::simulate(){
	printf("initialization\n");
	initializeProfile();

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
	fopen_s(&outShockWave, "./output/shock_wave.dat","w");
	double shockWaveR = 0;
	if(shockWavePoint >= 0 && shockWavePoint <= rgridNumber){
		shockWaveR = grid[shockWavePoint];
	}

	fprintf(outShockWave, "%d %lf %d %lf\n", 0, myTime, shockWavePoint, shockWaveR, 0, 0);
	fclose(outShockWave);
	deltaT = min2(minT, deltaT);

	currentIteration = 0;

	//основной цикл
	riemanSolve();

	while(myTime < maxTime && currentIteration < iterationNumber){
		++currentIteration;
		printf("iteration є %d\n", currentIteration);
		printf("time = %lf\n", myTime);
		printf("solving\n");

		evaluateHydrodynamic();
		riemanMove();

		myTime = myTime + deltaT;

		updateAll();

		updateMaxSoundSpeed();
		updateShockWavePoint();
		updateParameters();

		updateTimeStep();
		deltaT = min2(minT, deltaT);
		if(currentIteration % writeParameter == 0){
			//вывод на некоторых итераци€х
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

	//updateFluxes();

	TracPen(tempDensity, dFlux, 0);

	TracPen(tempMomentum, mFlux, 0);

	TracPen(tempEnergy, eFlux, 0);
}


//расчет разрывов, задача –имана

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


//“рак и ѕен, немного модифицированные

void Simulation::TracPen(double* u, double* flux, double cs){

	tempU[0] = u[0] - deltaT*(flux[1] - flux[0])/deltaR[0];
	for(int i = 1; i < rgridNumber - 1; ++i){
		tempU[i] = u[i] - deltaT*(flux[i+1] - flux[i])/deltaR[i];
		alertNaNOrInfinity(tempU[i],"tempU = NaN");
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

		double lambda1 = pointVelocity[i] - pointSoundSpeed[i];
		double lambda2 = pointVelocity[i];
		double lambda3 = pointVelocity[i] + pointSoundSpeed[i];

		//lambdaMod[0] = min2(middleVelocity[i-1] - sqrt(gamma*middlePressure[i-1]/middleDensity[i-1]), pointVelocity[i] - pointSoundSpeed[i]);
		lambdaMod[0] = lambda1;
		lambdaMod[1] = lambda2;
		lambdaMod[2] = lambda3;
		//lambdaMod[2] = max2(middleVelocity[i] + sqrt(gamma*middlePressure[i]/middleDensity[i]), pointVelocity[i] + pointSoundSpeed[i]);


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

	deltaT = 0.9*tempdt;
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
			tempShockWavePoint = i+1;
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
	for(int i = 0; i < rgridNumber; ++i){
		mass += middleDensity[i]*volume(i);
		totalMomentum += momentum(i)*volume(i);
		totalKineticEnergy += kineticEnergy(i)*volume(i);
		totalTermalEnergy += termalEnergy(i)*volume(i);
	}
	mass -= deltaT*(middleDensity[0]*middleVelocity[0] - middleDensity[rgridNumber-1]*middleVelocity[rgridNumber-1]);
	totalMomentum -= deltaT*(middleDensity[0]*sqr(middleVelocity[0]) + middlePressure[0] - middleDensity[rgridNumber-1]*sqr(middleVelocity[rgridNumber-1]) - middlePressure[rgridNumber-1]);

	totalKineticEnergy -= deltaT*(middleDensity[0]*cube(middleVelocity[0]) - middleDensity[rgridNumber-1]*cube(middleVelocity[rgridNumber-1]))/2;
	totalTermalEnergy -= deltaT*(middlePressure[0]*middleVelocity[0] - middlePressure[rgridNumber-1]*middleVelocity[rgridNumber-1])*gamma/(gamma-1);
	totalEnergy = totalTermalEnergy + totalKineticEnergy;
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
}

void Simulation::riemanSolve(){
	int i = shockWavePoint;
	double p;
	p1 = middlePressure[i-1];
	p6 = middlePressure[i];
	rho1 = middleDensity[i-1];
	rho6 = middleDensity[i];
	u1 = middleVelocity[i-1];
	u6 = middleVelocity[i];
	double c1 = sqrt(gamma*p1/rho1);
	double c6 = sqrt(gamma*p6/rho6);
	double u;
	double R3;
	double R4;
	double alpha1;
	double alpha6;

	double Uvacuum = -2*c1/(gamma - 1) - 2*c6/(gamma - 1);
		
	successiveApproximationPressure(p, u, R3, R4, alpha1, alpha6, p1, p6, u1, u6, rho1, rho6);

	bool isLeftShockWave = (p > p1);
	bool isRightShockWave = (p > p6);

	velocityCD = u;

	double D1;
	double D2;

	if(isLeftShockWave){
		D1 = u1 - alpha1/rho1;
	} else {
		D1 = u1 - c1;
	}

	if(isRightShockWave){
		D2 = u6 + alpha6/rho6;
	} else {
		D2 = u6 + c6;
	}

	velocityL1 = D1;
	velocityR1 = D2;

	velocityL2 = velocityL1;
	velocityR2 = velocityR1;

	rho3 = R3;
	u3 = u;
	p3 = p;
	rho4 = R4;
	u4 = u;
	p4 = p;

	p2 = p3;
	u2 = u3;
	rho2 = rho3;

	p5 = p4;
	u5 = u3;
	rho5 = rho4;

	if(!isLeftShockWave){
		velocityL2 = u - c1 - (gamma - 1)*(u1 - u)/2;
		u2 = (gamma - 1)*middleVelocity[i-1]/(gamma + 1) + 2*c1/(gamma + 1);
		p2 = middlePressure[i-1]*power(abs(u/c1), 2*gamma/(gamma - 1));
		rho2 = gamma*p/(u*u);
	} 
			
	
	if(!isRightShockWave){
		velocityR2 = u + c6 - (gamma - 1)*(u6 - u)/2;
		u5 = (gamma - 1)*middleVelocity[i]/(gamma + 1) - 2*c6/(gamma + 1);
		p5 = middlePressure[i]*power(abs(u/c6), 2*gamma/(gamma - 1));
		rho5 = gamma*p/(u*u);
	}

	if(velocityR1 < velocityR2){
		printf("incorrect velocity\n");
	}

	if(velocityR2 < velocityCD){
		printf("incorrect velocity\n");
	}

	if(velocityCD < velocityL2){
		printf("incorrect velocity\n");
	}

	if(velocityL2 < velocityL1){
		printf("incorrect velocity\n");
	}
}

void Simulation::riemanMove(){
	rCD += deltaT*velocityCD;
	rL1 += deltaT*velocityL1;
	rL2 += deltaT*velocityL2;
	rR1 += deltaT*velocityR1;
	rR2 += deltaT*velocityR2;
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