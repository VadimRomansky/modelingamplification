#include <list>
#include <time.h>
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

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
	double C = middlePressure[0]/power(middleDensity[0], _gamma);
	for(int i = 1; i < shockWavePoint; ++i){
		double momentumF = middleDensity[0]*middleVelocity[0]*middleVelocity[0] + middlePressure[0] - (cosmicRayPressure[i] - cosmicRayPressure[0]) - (magneticEnergy[i] - magneticEnergy[0])/2;

		double u = newton(i, momentumF);
		if(u < 0){
			printf("u < 0\n");
		}

		tempDensity[i] = middleDensity[0]*middleVelocity[0]/u;
		if(tempDensity[i] < 0){
			printf("density < 0\n");
			tempDensity[i] = epsilon*density0;
		}
		tempMomentum[i] = tempDensity[i]*u;
		tempEnergy[i] = tempDensity[i]*u*u*0.5 + C*power(tempDensity[i], _gamma)/(_gamma-1);
		if(tempEnergy[i] < 0){
			printf("energy < 0\n");
			tempEnergy[i] = tempDensity[i]*u*u/2;
		}

		alertNaNOrInfinity(tempDensity[i], "density = NaN");
		alertNaNOrInfinity(tempMomentum[i], "momentum = NaN");
		alertNaNOrInfinity(tempEnergy[i], "energy = NaN");
	}

	double p1 = middlePressure[shockWavePoint-1]/(middleDensity[0]*middleVelocity[0]*middleVelocity[0]);
	double pw1 = 0.5*magneticEnergy[shockWavePoint-1]/(middleDensity[0]*middleVelocity[0]*middleVelocity[0]);
	double u1 = tempMomentum[shockWavePoint-1]/(middleVelocity[0]*middleDensity[shockWavePoint-1]);
	double lambda = 2*_gamma*(p1 + pw1) + u1*(_gamma - 1);
	double Rtot = (lambda - sqrt2(lambda*lambda - 8*pw1*u1*(_gamma+1)*(_gamma-2)))/(4*u1*pw1*(_gamma-2));
	if(pw1 < 1E-10){
		Rtot = (_gamma + 1)/lambda;
	}
	double Rsub = Rtot*middleDensity[0]/middleDensity[shockWavePoint-1];

	if(Rsub < 1){
		printf("Rsub < 1\n");
	}

	double delta = sqr(Rsub - 1)*pw1/p1;
		
	double p2 = p1*((_gamma+1)*Rsub - (_gamma-1) + (_gamma-1)*(Rsub-1)*delta)/((_gamma+1) - (_gamma-1)*Rsub);

	if(abs(Rsub - (_gamma+1)/(_gamma-1)) < 1E-12){
		p2 = (middleDensity[shockWavePoint-1]*middleVelocity[shockWavePoint-1]*middleVelocity[shockWavePoint-1])*(Rsub-1)/Rsub;
	}

	alertNaNOrInfinity(p2, "p2 = NaN");

	if(p2 < 0){
		printf("p < 0\n");
	}

	/*double momentumF = middleDensity[0]*middleVelocity[0]*middleVelocity[0] + middlePressure[0] - (magneticEnergy[shockWavePoint] - magneticEnergy[0])/2;
	double energyF = (energy(shockWavePoint-1) + middlePressure[shockWavePoint-1])*middleVelocity[shockWavePoint-1] - 0.5*magneticEnergy[shockWavePoint]*(middleVelocity[shockWavePoint] - middleVelocity[shockWavePoint-1]);
	for(int k = 0; k < kgridNumber; ++k){
		energyF += growth_rate[shockWavePoint][k]*magneticField[shockWavePoint][k]*kgrid[k]*deltaLogK*middleDeltaR[shockWavePoint];
	}

	double a = middleDensity[0]*middleVelocity[0]*(0.5 - _gamma/(_gamma-1));
	double b = _gamma*momentumF/(_gamma - 1);
	double c = - energyF;

	double D = b*b - 4*a*c;
	if(D < 0){
		printf("discriminant < 0 i = %d\n",shockWavePoint);
		D = 0;
	}
	double u2 = (-b - sqrt(D))/(2*a);
	double u1 = (-b + sqrt(D))/(2*a);

	double d = tempMomentum[shockWavePoint-1]/u1;
	if(d < 0){
		printf("density < 0\n");
		d = epsilon*density0;
	}
	double m = tempMomentum[shockWavePoint-1];
	double e = d*u1*u1*(0.5 - 1/(_gamma-1)) + momentumF/(_gamma - 1);
	if(e < 0){
		printf("energy < 0\n");
		e = d*u1*u1/2;
	}*/

	double d = middleDensity[0]*Rtot;
	double m = middleDensity[0]*middleVelocity[0];
	double u = m/d;
	double e = d*u*u/2 + p2/(_gamma-1);


	for(int i = shockWavePoint; i < rgridNumber; ++i){
		tempDensity[i] = d;
		tempMomentum[i] = m;
		tempEnergy[i] = e;

		alertNaNOrInfinity(tempDensity[i], "density = NaN");
		alertNaNOrInfinity(tempMomentum[i], "momentum = NaN");
		alertNaNOrInfinity(tempEnergy[i], "energy = NaN");
	}
}

double Simulation::newton(int i, double momentumF){
	double u = middleVelocity[i];
	double startU = 0.8*u;
	double endU = 1.2*u;
	double C = middlePressure[0]/power(middleDensity[0], _gamma);
	double rho0 = middleDensity[0];
	double u0 = middleVelocity[0];
	for(int i = 0; i < 30; ++i){
		double phi = u - (momentumF/(rho0*u0)) + C*power(rho0*u0, _gamma-1)/power(u, _gamma);
		if(phi < 0){
			startU = u;
		} else {
			endU = u;
		}
		u = (startU + endU)/2;
		if(u < 0){
			printf("u < 0 in newton\n");
		}
	}
	/*for(int i = 0; i < 10; ++i){
		double phi = (momentumF/(rho0*u0)) - C*power(rho0*u0, _gamma-1)/power(u, _gamma);
		double dphi = C*_gamma*power(rho0*u0, _gamma-1)/power(u, _gamma+1);
		u = u - phi/dphi;
		if(u < 0){
			printf("u < 0 in newton\n");
		}
	}*/
	return u;
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
	for(int i = 0; i < rgridNumber; ++i){
		double x = 0.97;
		tempDensity[i] = x*middleDensity[i] + (1 - x)*tempDensity[i];
		tempMomentum[i] = x*momentum(i) + (1 - x)*tempMomentum[i];
		tempEnergy[i] = x*energy(i) + (1 - x)*tempEnergy[i];

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