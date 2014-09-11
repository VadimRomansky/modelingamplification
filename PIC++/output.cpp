#include "stdio.h"
#include "math.h"
#include "vector"

#include "output.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"

void outputDistribution(FILE* outFile, std::vector<Particle*> particles, int particleType){
	double minMomentum = -1; //todo something else
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
	fprintf(outFile, "%g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n", time, particle->coordinates.x, particle->coordinates.y, particle->coordinates.z, particle->momentum.x, particle->momentum.y, particle->momentum.z, particle->momentum.getNorm());
}

void outputGrid(FILE* outFile, double* grid, int number) {
	for(int i = 0; i < number; ++i) {
		fprintf(outFile, "%15.10g %15.10g\n", grid[i], (grid[i+1] + grid[i])/2);
	}
}

void outputFields(FILE* outEfile, FILE* outBfile, Vector3d*** Efield, Vector3d*** Bfield, int xnumber, int ynumber, int znumber) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fprintf(outEfile, "%15.10g %15.10g %15.10g\n", Efield[i][j][k].x, Efield[i][j][k].y, Efield[i][j][k].z);
				fprintf(outBfile, "%15.10g %15.10g %15.10g\n", Bfield[i][j][k].x, Bfield[i][j][k].y, Bfield[i][j][k].z);
			}
		}
	}
}