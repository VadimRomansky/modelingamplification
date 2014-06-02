#ifndef INPUT_H
#define INPUT_H
#include "simulation.h"
#include "constants.h"

Simulation* readInput();
void deserialize(FILE* hydroFile, FILE* distributionFile, FILE* fieldFile, FILE* gridFile, FILE* pgridFile, FILE* kgridFile, FILE* infoFile, Simulation* simulation);

#endif