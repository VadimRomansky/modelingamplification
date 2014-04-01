#include <time.h>
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"

double Simulation::diffusionCoef(double x, double p){
	double B = B0;
	return p*speed_of_light*speed_of_light/(electron_charge*B);
}

//инжекционный член
double Simulation::injection(){
	double pf = pgrid[injectionMomentum];
	double dp = (pgrid[injectionMomentum + 1] - pgrid[injectionMomentum - 1])/2;
	return 0.1*density(0)*abs(velocity(0))/(pf*pf*massProton*dp);
	//return 1E-30;
}

//расчет космических лучей

void Simulation::evaluateCR(){

	printf("solve CR\n");
	/*if(shockWavePoint > 0 && shockWavePoint < rgridNumber){
		distributionFunction[shockWavePoint][injectionMomentum] += injection()*deltaT;
	}*/
	double* upper = new double[rgridNumber];
	double* middle = new double[rgridNumber];
	double* lower = new double[rgridNumber];

	double* f = new double[rgridNumber];
 	double* solution = new double[rgridNumber];

	double* alpha = new double[rgridNumber];
	double* beta = new double[rgridNumber];

	for(int j = 0; j < pgridNumber; ++j){

		double y = logPgrid[j];
		double p = pgrid[j];
		double functionMomentumRight = distributionFunction[0][j];
		double functionMomentumLeft = 0;
		if(j > 0){
			functionMomentumLeft = distributionFunction[0][j-1];
		}

		double dx = (grid[1] - grid[0]);
		double dxleft = grid[1] - grid[0];
		double dxright = grid[1] - grid[0];
		double xleft = grid[0] - dxleft/2;
		double xright = (grid[1] + grid[0])/2;

		double functionXRight = (distributionFunction[1][j] + distributionFunction[0][j])/2;
		double functionXLeft = (distributionFunction[0][j]/2);

		middle[0] = 1 + (deltaT/(2*dx))*(diffusionCoef(xright,p)/dxright + diffusionCoef(xleft,p)/dxleft); 
		f[0] = distributionFunction[0][j] + (deltaT/(2*dx))*(diffusionCoef(xright,p)*distributionFunction[1][j]/dxright) - (deltaT/(2*dx))*(diffusionCoef(xleft, p)/dxleft + diffusionCoef(xright, p)/dxright)*distributionFunction[0][j]
		-(deltaT/dx)*(velocity(xright)*functionXRight - velocity(xleft)*functionXLeft) + (deltaT/3)*((velocity(xright) - velocity(xleft))/dx)*(functionMomentumRight -functionMomentumLeft)/deltaLogP;
		lower[0] = 0;
		upper[0] = -(deltaT/(2*dx))*diffusionCoef(xright,p)/dxright;
		upper[rgridNumber - 1] = 0;

		for(int i = 1; i < rgridNumber; ++i){
			functionMomentumRight = distributionFunction[i][j];
			functionMomentumLeft = 0;
			if(j > 0){
				functionMomentumLeft = distributionFunction[i][j-1];
			}
			dx = (grid[i+1] - grid[i-1])/2;
			dxright = grid[i+1] - grid[i];
			dxleft = grid[i] - grid[i-1];
			xright = (grid[i+1] + grid[i])/2;
			xleft = (grid[i] + grid[i-1])/2;

			functionXRight = (distributionFunction[i][j] + distributionFunction[i+1][j])/2;
			functionXLeft = (distributionFunction[i][j] + distributionFunction[i-1][j])/2;

			lower[i] = -(deltaT/(2*dx))*(diffusionCoef(xleft, p)/dxleft);
			alertNaNOrInfinity(lower[i],"lower = NaN");
			upper[i] = -(deltaT/(2*dx))*(diffusionCoef(xright, p)/dxright);
			alertNaNOrInfinity(upper[i],"lower = NaN");
			middle[i] = 1 - lower[i] - upper[i];
			f[i] = distributionFunction[i][j] - upper[i]*(distributionFunction[i+1][j] - distributionFunction[i][j])
											  + lower[i]*(distributionFunction[i][j] - distributionFunction[i-1][j])
											  - (deltaT/dx)*(velocity(xright)*functionXRight - velocity(xleft)*functionXLeft)
											  //- (deltaT/dx)*(velocity(xright)*distributionFunction[i][j] - velocity(xleft)*distributionFunction[i-1][j])
											  + (deltaT/3)*((velocity(xright) - velocity(xleft))/dx)*((functionMomentumRight - functionMomentumLeft)/deltaLogP);
			alertNaNOrInfinity(f[i],"f[i] = NaN");
			alertNegative(f[i],"f[i] < 0");
		}
		if(j == injectionMomentum){
			f[rgridNumber/2] += deltaT*injection();
		}
		
		solveThreeDiagonal(middle, upper, lower, f, solution, alpha, beta);
		
		for(int i = 0; i < rgridNumber; ++i){
			tempDistributionFunction[i][j] = solution[i];
			if(solution[i] < 0){
				if(abs(solution[i]) > 1E-300){
					tempDistributionFunction[i][j] = 0;
					printf("tempDistribution < 0\n");
				} else {
					tempDistributionFunction[i][j] = 0;
				}
			}
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
	delete[] solution;
	delete[] alpha;
	delete[] beta;
}

//решение трёх диагональной матрицы
void Simulation::solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x, double* alpha, double* beta){

	alpha[1] = -upper[0]/middle[0];
	beta[1] = f[0]/middle[0];
	for(int i = 2; i < rgridNumber; ++i){
		double temp = lower[i-1]*alpha[i-1] + middle[i-1];
		alpha[i] = -upper[i-1]/temp;
		beta[i] = (f[i-1] - lower[i-1]*beta[i-1])/temp;
	}

	x[rgridNumber - 1] = (f[rgridNumber-1] - lower[rgridNumber-1]*beta[rgridNumber-1])/(lower[rgridNumber-1]*alpha[rgridNumber-1] + middle[rgridNumber-1]);
	alertNaNOrInfinity(x[rgridNumber-1],"x = NaN");
	//alertNegative(x[rgridNumber-1],"x < 0");

	for(int i = rgridNumber - 2; i >= 0; --i){
		x[i] = alpha[i+1]*x[i+1] + beta[i+1];
		alertNaNOrInfinity(x[i],"x = NaN");
		//alertNegative(x[i],"x < 0");
	}
}