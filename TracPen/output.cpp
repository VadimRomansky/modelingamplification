#include "stdafx.h"
#include "output.h"
#include "constants.h"

void output(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		fprintf(outFile,"%17.12lf %17.12lf %38.30lf %17.12lf %17.12lf\n", simulation->bins[i]->r, simulation->bins[i]->U, simulation->bins[i]->density, simulation->bins[i]->pressure, simulation->bins[i]->temperature);
	}
}