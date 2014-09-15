#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include "vector"
#include "stdio.h"

#include "vector3d.h"
#include "particle.h"

void outputDistribution(FILE* outFile, std::vector<Particle*> particles, int particleType);
void outputTraectory(FILE* outFile, Particle* particle, double time);
void outputGrid(FILE* outFile, double* grid, int number);
void outputFields(FILE* outEfile, FILE* outBfile, Vector3d*** Efield, Vector3d*** Bfield, int xnumber, int ynumber, int znumber, double plasma_priod, double gyroradius);


#endif