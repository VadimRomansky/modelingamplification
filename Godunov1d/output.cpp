#include "stdafx.h"
#include "output.h"
#include "constants.h"

void output(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		fprintf(outFile,"%17.12lf %17.12lf %38.30lf %28.20lf %17.12lf\n", simulation->grid[i], simulation->pointVelocity[i], simulation->pointDensity[i], simulation->pointPressure[i], simulation->temperatureAtPoint(i));
	}
}