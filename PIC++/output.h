#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include "vector"
#include "stdio.h"

#include "vector3d.h"
#include "particle.h"
#include "simulation.h"

void outputDistribution(FILE* outFile, std::vector<Particle*> particles, int particleType);
void outputTraectory(FILE* outFile, Particle* particle, double time);
void outputGrid(FILE* outFile, double* grid, int number);
void outputFields(FILE* outEfile, FILE* outBfile, double*** EfieldX, double*** EfieldY, double*** EfieldZ, double*** BfieldX, double*** BfieldY, double*** BfieldZ, int xnumber, int ynumber, int znumber, double plasma_priod, double gyroradius);
void outputConcentrations(FILE* outFile, double*** electronConcentration, double*** protonConcentration, double*** chargeDensity, int xnumber, int ynumber, int znumber, double plasma_period, double gyroradius);
void outputVelocity(FILE* outFile, Vector3d*** velocity, int xnumber, int ynumber, int znumber, double plasma_period, double gyroradius);
void outputArrayParameter(FILE* outFile, double*** arrayParameter, int xnumber, int ynumber, int znumber);
void outputGeneral(FILE* outFile, Simulation* simulatiom);
void outputDivergenceError(FILE* outFile, Simulation* simulation);


#endif