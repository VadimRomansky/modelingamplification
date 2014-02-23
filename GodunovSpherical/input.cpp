#include <stdio.h>
#include "stdafx.h"
#include "input.h"
#include "simulation.h"

Simulation* readInput(){
	FILE* infile = fopen("tamc_in.dat","r");
	if (infile != NULL){
		char* s = NULL;
		Simulation* simulation = new Simulation();
		fscanf(infile,"%d",&simulation->iterationNumber);
        fscanf(infile,"%d",&simulation->particlesNumber);
		fscanf(infile,"%lf",&simulation->U0);
		fscanf(infile,"%lf",&simulation->density0);
        fscanf(infile,"%lf",&simulation->B0);
        fscanf(infile,"%lf",&simulation->temperature);
		fscanf(infile,"%lf",&simulation->upstreamR);
		fscanf(infile,"%d",&simulation->rgridNumber);
		fscanf(infile,"%d",&simulation->simulationType);
		fscanf(infile,"%lf",&simulation->maxTime);
		fclose(infile);
		return simulation;
	} else {
		return NULL;
	}
}
	