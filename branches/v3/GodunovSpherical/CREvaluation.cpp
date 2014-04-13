#include <time.h>
#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "output.h"
#include "progon.h"

double Simulation::diffusionCoef(int i, double p){
	double B = B0;
	double coef = 10000*p*speed_of_light*speed_of_light/(electron_charge*B);
	double dx = deltaR[i];
	double lambda = coef/speed_of_light;
	return coef;
}

//������������ ����
double Simulation::injection(){
	double pf = pgrid[injectionMomentum];
	double dp = (pgrid[injectionMomentum + 1] + pgrid[injectionMomentum - 1])/2;
	return middleDensity[shockWavePoint]*abs(middleVelocity[shockWavePoint])/(pf*pf*massProton*dp);
	//return 1;
}


//������ ����������� �����

void Simulation::evaluateCR(){
	printf("solve CR\n");
	/*if(shockWavePoint > 0 && shockWavePoint < rgridNumber){
		distributionFunction[shockWavePoint][injectionMomentum] += injection()*deltaT;
	}*/
	
	double* upper = new double[rgridNumber+1];
	double* middle = new double[rgridNumber+1];
	double* lower = new double[rgridNumber+1];

	double* f = new double[rgridNumber+1];
 	double* x = new double[rgridNumber+1];

	double* alpha = new double[rgridNumber];
	double* beta = new double[rgridNumber];

	for(int k = 0; k < pgridNumber; ++k){
		double y = logPgrid[k];
		double p = pgrid[k];
		double gkp = distributionFunction[0][k];
		double gkm=0.0;
		if (k == 0) {
			gkm=0.0;
		} else {
			gkm = distributionFunction[0][k-1];
		}
		double gkpp;
		if(k < pgridNumber-1){
			gkpp = distributionFunction[0][k+1];
		} else {
			gkpp = 0;
		}
		double dx = (grid[1] + upstreamR/2)/2;
		double dxp = grid[1] - grid[0];
		double dxm = grid[0] + upstreamR/2;
		double xp = (grid[1] + grid[0])/2;
		double xm = (grid[0] - upstreamR/2)/2;
		double gip = (distributionFunction[1][k] + distributionFunction[0][k])/2;
		double gim = distributionFunction[0][k]/2;
		double dV;
		middle[0] = 1;// + (deltaT/(2*dx))*(diffusionCoef(0,p)/dxp + diffusionCoef(0,p)/dxm);
		upper[0] = 0;// -(deltaT/(2*dx))*(diffusionCoef(0,p)/dxp);
		f[0]=distributionFunction[0][k];// + (deltaT/(2*dx))*(diffusionCoef(0,p)*distributionFunction[1][k]/dxp)
					 // - (deltaT/(2*dx))*(diffusionCoef(0,p)/dxp+diffusionCoef(0,p)/dxm)*distributionFunction[0][k] 
					 // - (deltaT/dx)*(middleVelocity[1]*gip - middleVelocity[0]*gim)
					 // + (deltaT/3)*((middleVelocity[1] - middleVelocity[0])/dx)*((gkp-gkm)/deltaLogP);
		for(int i = 1; i < rgridNumber-1; ++i){
			gkp = distributionFunction[i][k];
			if (k == 0) {
				gkm=0;
			} else {
				gkm=distributionFunction[i][k-1];
			}
			if(k < pgridNumber-1){
				gkpp = distributionFunction[i][k+1];
			} else {
				gkpp = 0;
			}
			dx = (grid[i+1] - grid[i-1])/2;
			dxp=grid[i+1]-grid[i];
			dxm=grid[i]-grid[i-1];
			xp=(grid[i+1]+grid[i])/2;
			xm=(grid[i]+grid[i-1])/2;
			dV = (xp*xp*xp - xm*xm*xm)/3;
			gip=(distributionFunction[i+1][k] + distributionFunction[i][k])/2;
			gim=(distributionFunction[i][k] + distributionFunction[i-1][k])/2;
			lower[i-1] = -(deltaT/(2*dV))*(xm*xm*diffusionCoef(i-1,p)/dxm);
			middle[i] = 1 + (deltaT/(2*dV))*(xp*xp*diffusionCoef(i-1,p)/dxp + xm*xm*diffusionCoef(i,p)/dxm);
			upper[i] = -(deltaT/(2*dV))*(xp*xp*diffusionCoef(i,p)/dxp);
			f[i] = distributionFunction[i][k] + (deltaT/(2*dV))*(xp*xp*diffusionCoef(rgridNumber-1,p)*(distributionFunction[i+1][k] - distributionFunction[i][k])/dxp
							- xm*xm*diffusionCoef(i-1,p)*(distributionFunction[i][k] - distributionFunction[i-1][k])/dxm)
							- (deltaT/dV)*(xp*xp*middleVelocity[i]*distributionFunction[i][k] - xm*xm*middleVelocity[i-1]*distributionFunction[i-1][k]);
			if((xp*xp*middleVelocity[i] - xm*xm*middleVelocity[i-1]) < 0){
				if(gkp - gkm < 0)
					f[i] += (deltaT/3)*((xp*xp*middleVelocity[i] - xm*xm*middleVelocity[i-1])/dV)*((gkp - gkm)/deltaLogP);
			} else {
				if(gkpp - gkp > 0)
					f[i] += (deltaT/3)*((xp*xp*middleVelocity[i] - xm*xm*middleVelocity[i-1])/dV)*((gkpp - gkp)/deltaLogP);
			}
			if(i == shockWavePoint && k == injectionMomentum){
				f[i] += deltaT*injection()*grid[i]*grid[i]/dV;
				injectedParticles += injection()*deltaT*4*pi*volume(i)*deltaLogP*grid[i]*grid[i]/dV;
			}
		}
		gkp = distributionFunction[rgridNumber-1][k];
		if (k==0) {
			gkm=0;
		} else {
			gkm = distributionFunction[rgridNumber-1][k-1];
		}
		dx = (grid[rgridNumber-1] - grid[rgridNumber - 2])/2;
		dxm = (grid[rgridNumber-1] - grid[rgridNumber - 2]);
		xp = grid[rgridNumber-1];
		xm = (grid[rgridNumber-1] + grid[rgridNumber-2])/2;
		dV = (xp*xp*xp - xm*xm*xm)/3;
		gip = distributionFunction[rgridNumber-1][k];
		gim = (distributionFunction[rgridNumber-1][k] + distributionFunction[rgridNumber - 2][k])/2;
		lower[rgridNumber-2] = -(deltaT/(2*dV))*(xp*xp*diffusionCoef(rgridNumber-1,p)/dxp);
		middle[rgridNumber-1] = 1 + (deltaT/(2*dV))*(xm*xm*diffusionCoef(rgridNumber-1,p)/dxm);
		f[rgridNumber-1] = distributionFunction[rgridNumber-1][k] - (deltaT/(2*dV))*(xm*xm*diffusionCoef(rgridNumber-2,p)*(distributionFunction[rgridNumber-1][k] - distributionFunction[rgridNumber - 2][k])/dxm)
							  - (deltaT/dV)*grid[rgridNumber-1]*grid[rgridNumber-1]*(middleVelocity[rgridNumber-1]*gip - middleVelocity[rgridNumber-2]*gim)
							  + (deltaT/3)*((xp*xp*middleVelocity[rgridNumber-1] - xm*xm*middleVelocity[rgridNumber-2])/dV)*((gkp - gkm)/deltaLogP);
		progon(lower,middle, upper,rgridNumber-1,f,x, alpha, beta);

		for(int i = 0; i < rgridNumber; ++i){
			//alertNegative(x[i],"tempDistribution < 0");
			alertNaNOrInfinity(x[i],"tempDistribution <= NaN");
			tempDistributionFunction[i][k]= x[i];
			if(x[i] < 0){
				tempDistributionFunction[i][k] = 0;
				if(abs(x[i]) > 1E-50){
					printf("distribution[i] < 0\n");
				}
			}
		}
		//tempDistributionFunction[rgridNumber][k] =x[rgridNumber-1];
		tempDistributionFunction[rgridNumber][k] = 0;
	}

	for(int i = 0; i <= rgridNumber; ++i){
		for(int j = 0; j < pgridNumber; ++j){
			distributionFunction[i][j] = tempDistributionFunction[i][j];
		}
	}
	evaluateCosmicRayPressure();

	delete[] upper;
	delete[] middle;
	delete[] lower;
	delete[] f;
	delete[] x;
	delete[] alpha;
	delete[] beta;
}

//������� ��� ������������ �������
void Simulation::solveThreeDiagonal(double* middle, double* upper, double* lower, double* f, double* x, double* alpha, double* beta){

	alpha[1] = -upper[0]/middle[0];
	beta[1] = f[0]/middle[0];
	for(int i = 2; i < rgridNumber; ++i){
		double temp = lower[i-1]*alpha[i-1] + middle[i-1];
		alpha[i] = -upper[i-1]/temp;
		beta[i] = (f[i-1] - lower[i-1]*beta[i-1])/temp;
	}

	x[rgridNumber - 1] = (f[rgridNumber-1] - lower[rgridNumber-1]*beta[rgridNumber-1])/(lower[rgridNumber-1]*alpha[rgridNumber-1] + middle[rgridNumber-1]);
	//alertNaNOrInfinity(x[rgridNumber-1],"x = NaN");
	//alertNegative(x[rgridNumber-1],"x < 0");

	for(int i = rgridNumber - 2; i >= 0; --i){
		x[i] = alpha[i+1]*x[i+1] + beta[i+1];
		//alertNaNOrInfinity(x[i],"x = NaN");
		//alertNegative(x[i],"x < 0");
	}
}


//���������� �������� ����������� �����

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
		pressure *= 4*pi;
		cosmicRayPressure[i] = pressure;
	}

	delete[] partPressure;
}