#include <stdio.h>
#include "stdafx.h"
#include "input.h"
#include "simulation.h"

Simulation* readInput(){
	FILE* infile = fopen("tamc_in.dat","r");
	if (infile != NULL){
		std::string s;
		Simulation* simulation = new Simulation();
		fscanf(infile,"%d",&simulation->iterationNumber);
        fscanf(infile,"%d%",&simulation->particlesNumber);
		fscanf(infile,"%lf",&simulation->U0);
		fscanf(infile,"%lf",&simulation->Rtot);
		fscanf(infile,"%lf",&simulation->density0);
        fscanf(infile,"%lf",&simulation->B0);
        fscanf(infile,"%lf",&simulation->temperature);
		fscanf(infile,"%lf",&simulation->X1);
		fscanf(infile,"%lf",&simulation->X2);
		fscanf(infile,"%lf",&simulation->Y1);
		fscanf(infile,"%lf",&simulation->Y2);
		int a;
		fscanf(infile,"%d", &a);
		simulation->smallAngleScattering = (a == 1);
		fscanf(infile,"%d",&simulation->xgridNumber);
		fscanf(infile,"%d",&simulation->ygridNumber);
		fscanf(infile,"%d",&simulation->shockWaveIndex);
		fscanf(infile,"%d",&simulation->freeTimeEvaluationType);
		fscanf(infile,"%d",&simulation->averageVelocityEvaluationType);
		fclose(infile);
		return simulation;
	} else {
		return NULL;
	}
}
	