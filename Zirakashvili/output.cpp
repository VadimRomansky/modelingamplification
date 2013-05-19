#include "stdafx.h"
#include "output.h"
#include "constants.h"

void output(FILE* outFile, Simulation* simulation){
	for(int i = 0; i < simulation->rgridNumber; ++i){
		fprintf(outFile,"%lf %lf %lf %lf %lf\n", simulation->upstreamBins1[i]->r, simulation->upstreamBins1[i]->U, simulation->upstreamBins1[i]->density, simulation->upstreamBins1[i]->pressure, simulation->upstreamBins1[i]->temperature);
	}
	for(int i = 0; i < simulation->rgridNumber; ++i){
		fprintf(outFile,"%lf %lf %lf %lf %lf\n", simulation->downstreamBins1[i]->r, simulation->downstreamBins1[i]->U, simulation->downstreamBins1[i]->density, simulation->downstreamBins1[i]->pressure, simulation->downstreamBins1[i]->temperature);
	}
	for(int i = 0; i < simulation->rgridNumber; ++i){
		fprintf(outFile,"%lf %lf %lf %lf %lf\n", simulation->downstreamBins2[i]->r, simulation->downstreamBins2[i]->U, simulation->downstreamBins2[i]->density, simulation->downstreamBins2[i]->pressure, simulation->downstreamBins2[i]->temperature);
	}
	for(int i = 0; i < simulation->rgridNumber; ++i){
		fprintf(outFile,"%lf %lf %lf %lf %lf\n", simulation->upstreamBins2[i]->r, simulation->upstreamBins2[i]->U, simulation->upstreamBins2[i]->density, simulation->upstreamBins2[i]->pressure, simulation->upstreamBins2[i]->temperature);
	}
}