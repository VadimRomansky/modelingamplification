#include "stdio.h"
#include "math.h"
#include "vector"

#include "output.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"

void outputDistribution(FILE* outFile, std::vector<Particle*> particles, int particleType){
	double minMomentum = -1; //todo somethingelse
	double maxMomentum = 0;
	for(int i = 0; i < particles.size(); ++i){
		if(particles[i]->type == particleType){
			if(minMomentum < 0){
				minMomentum = particles[i]->momentum.getNorm();
			}
			if(particles[i]->momentum.getNorm() < minMomentum){
				minMomentum = particles[i]->momentum.getNorm();
			} else {
				if(particles[i]->momentum.getNorm() > maxMomentum){
					maxMomentum = particles[i]->momentum.getNorm();
				}
			}
		}
	}

	double pgrid[pnumber+1];
	double distribution[pnumber];
	double logMinMomentum = log(minMomentum);
	pgrid[0] = minMomentum;
	distribution[0] = 0;
	double deltaLogP = (log(maxMomentum) - log(minMomentum))/(pnumber);
	for(int i = 1; i < pnumber; ++i){
		distribution[i] = 0;
		pgrid[i] = exp(logMinMomentum + i*deltaLogP);
	}
	pgrid[pnumber] = maxMomentum;

	double weight = 0;

	for(int i = 0; i < particles.size(); ++i){
		if(particles[i]->type == particleType){
			int j = (log(particles[i]->momentum.getNorm()) - logMinMomentum)/deltaLogP;
			if( j >= 0 && j < pnumber){
				distribution[j] += particles[i]->weight;
				weight += particles[i]->weight;
			}
		}
	}

	for(int i = 0; i < pnumber; ++i){
		distribution[i] /= (weight*(pgrid[i+1] - pgrid[i]));
	}

	for(int i = 0; i < pnumber; ++i){
		fprintf(outFile, "%g %g\n", (pgrid[i] + pgrid[i+1])/2, distribution[i]);
	}
}

void outputTraectory(FILE* outFile, Particle* particle, double time){
	Vector3d velocity = particle->velocity();
	fprintf(outFile, "%g %g %g %g %g %g %g", time, particle->coordinates.x, particle->coordinates.y, particle->coordinates.z, velocity.x, velocity.y, velocity.z);
}