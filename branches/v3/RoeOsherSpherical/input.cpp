#include <stdio.h>
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

void deserialize(FILE* hydroFile, FILE* distributionFile, FILE* fieldFile, FILE* gridFile, FILE* pgridFile, FILE* kgridFile, FILE* infoFile, Simulation* simulation){
	fscanf(infoFile,"%d", &simulation->iterationNumber);
    fscanf(infoFile,"%d", &simulation->particlesNumber);
	fscanf(infoFile,"%lf", &simulation->U0);
	fscanf(infoFile,"%lf", &simulation->density0);
    fscanf(infoFile,"%lf", &simulation->B0);
    fscanf(infoFile,"%lf", &simulation->temperature);
	fscanf(infoFile,"%lf", &simulation->upstreamR);
	fscanf(infoFile,"%d", &simulation->rgridNumber);
	fscanf(infoFile,"%d", &simulation->simulationType);
	fscanf(infoFile,"%lf", &simulation->maxTime);
	fscanf(infoFile, "%lf", &simulation->myTime);
	fscanf(infoFile, "%d", &simulation->currentIteration);
	fscanf(infoFile, "%lf", &simulation->injectedParticles);

	for(int j = 0; j < pgridNumber; ++j){
		fprintf(pgridFile, "%lf", &simulation->pgrid[j]);
	}

	for(int j = 0; j < kgridNumber; ++j){
		fprintf(kgridFile, "%lf", &simulation->kgrid[j]);
	}

	for(int i = 0; i <= simulation->rgridNumber; ++i){
		fscanf(gridFile,"%lf", simulation->grid[i]);
		for(int j = 0; j < pgridNumber; ++j){
			fscanf(distributionFile, "%lf", &simulation->distributionFunction[i][j]);
		}
	}

	for(int i = 0; i < simulation->rgridNumber; ++i){
		fscanf(hydroFile, "%g %g %g", &simulation->middleDensity[i], &simulation->middleVelocity[i], &simulation->middlePressure[i]);
		for(int k = 0; k < kgridNumber; ++k){
			fscanf(fieldFile, "%g", &simulation->magneticField[i][k]);
		}
	}
}
	