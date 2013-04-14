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
		fscanf(infile,"%lf",&simulation->Rtot);
		fscanf(infile,"%lf",&simulation->density0);
        fscanf(infile,"%lf",&simulation->B0);
        fscanf(infile,"%lf",&simulation->temperature0);
		//fscanf(infile,"%lf",&simulation->freeEscapeBoundary);
		fscanf(infile,"%lf",&simulation->upstreamR);
		fscanf(infile,"%lf",&simulation->downstreamR);
		int a;
		fscanf(infile,"%d", &a);
		simulation->smallAngleScattering = (a == 1);
		fscanf(infile,"%d",&simulation->rgridNumber);
		fscanf(infile,"%d",&simulation->shockWavePoint);
		fclose(infile);
		return simulation;
	} else {
		return NULL;
	}
}
	