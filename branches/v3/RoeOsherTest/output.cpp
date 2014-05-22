#include "output.h"
#include "constants.h"
#include "util.h"

void output(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		double density;
		double velocity;
		double pressure;
		double r = simulation->middleGrid[i];

		if(r <= simulation->rL1){
			density = simulation->rho1;
			velocity = simulation->u1;
			pressure = simulation->p1;
		} else if(r <= simulation->rL2){
			density = simulation->rho2;
			velocity = simulation->u2;
			pressure = simulation->p2;
		} else if(r <= simulation->rCD){
			density = simulation->rho3;
			velocity = simulation->u3;
			pressure = simulation->p3;
		} else if(r <= simulation->rR2){
			density = simulation->rho4;
			velocity = simulation->u4;
			pressure = simulation->p4;
		} else if(r <= simulation->rR1){
			density = simulation->rho5;
			velocity = simulation->u5;
			pressure = simulation->p5;
		} else {
			density = simulation->rho6;
			velocity = simulation->u6;
			pressure = simulation->p6;
		}
		fprintf(outFile,"%g %g %g %g %g %g %g\n", simulation->middleGrid[i], simulation->middleVelocity[i], simulation->middleDensity[i], simulation->middlePressure[i], velocity, density, pressure);
	}
}


void outputNewGrid(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		fprintf(outFile, "%d %17.12lf\n", i, simulation->tempGrid[i]);
	}
}