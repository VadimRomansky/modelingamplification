#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include "vector"
#include "stdio.h"

#include "particle.h"

void outputDistribution(FILE* outFile, std::vector<Particle*> particles, int particleType);
void outputTraectory(FILE* outFile, Particle* particle, double time);

#endif