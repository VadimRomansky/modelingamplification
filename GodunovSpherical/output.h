#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include "simulation.h"

void output(FILE* outFile, Simulation* simulation);
void outputDistribution(FILE* distributionFile, FILE* fullDistributionFile, FILE* coordinateDistributionFile, Simulation* simulation);
void outputDistributionP3(FILE* distributionFile, FILE* fullDistributionFile, FILE* coordinateDistributionFile, Simulation* simulation);
void outputNewGrid(FILE* outFile, Simulation* simulation);
void outMatrix(double* a, double* c, double* b, int N, double* f, double* x);

#endif