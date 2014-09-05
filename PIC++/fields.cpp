#include "stdlib.h"
#include "stdio.h"

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"

void Simulation::evaluateFields(){
	updateParameters();

	evaluateMaxwellEquationMatrix();

	generalizedMinimalResidualMethod();

	evaluateMagneticField();

	for(int i = 0; i < xnumber + 1; ++i){
		for(int j = 0; j < ynumber + 1; ++j){
			for(int k = 0; k < znumber + 1; ++k){
				newEfield[i][j][k] = (tempEfield[i][j][k] - tempEfield[i][j][k]*(1 - theta))/theta;
			}
		}
	}
}

void Simulation::evaluateMaxwellEquationMatrix(){
}

void Simulation::evaluateMagneticField(){
}