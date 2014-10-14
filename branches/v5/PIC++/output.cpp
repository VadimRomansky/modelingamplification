#include "stdio.h"
#include "math.h"
#include "vector"

#include "output.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "util.h"

void outputDistribution(FILE* outFile, std::vector<Particle*> particles, int particleType){
	double minMomentum = -1; //todo something else
	double maxMomentum = 0;
	for(int i = 0; i < particles.size(); ++i){
		if(particles[i]->type == particleType){
			if(minMomentum < 0){
				minMomentum = particles[i]->momentum.norm();
			}
			if(particles[i]->momentum.norm() < minMomentum){
				minMomentum = particles[i]->momentum.norm();
			} else {
				if(particles[i]->momentum.norm() > maxMomentum){
					maxMomentum = particles[i]->momentum.norm();
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
			int j = (log(particles[i]->momentum.norm()) - logMinMomentum)/deltaLogP;
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
	fprintf(outFile, "%g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g\n", time, particle->coordinates.x, particle->coordinates.y, particle->coordinates.z, particle->momentum.x, particle->momentum.y, particle->momentum.z, particle->momentum.norm());
}

void outputGrid(FILE* outFile, double* grid, int number) {
	for(int i = 0; i < number; ++i) {
		fprintf(outFile, "%15.10g %15.10g\n", grid[i], (grid[i+1] + grid[i])/2);
	}
}

void outputFields(FILE* outEfile, FILE* outBfile, double*** EfieldX, double*** EfieldY, double*** EfieldZ, double*** BfieldX, double*** BfieldY, double*** BfieldZ, int xnumber, int ynumber, int znumber, double plasma_preiod, double gyroradius) {
	double scale = 1.0/(plasma_preiod*gyroradius);
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fprintf(outBfile, "%15.10g %15.10g %15.10g\n", scale*BfieldX[i][j][k], scale*BfieldY[i][j][k], scale*BfieldZ[i][j][k]);
			}
		}
	}

	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fprintf(outEfile, "%15.10g %15.10g %15.10g\n", scale*EfieldX[i][j][k], scale*EfieldY[i][j][k], scale*EfieldZ[i][j][k]);
			}
		}
	}
}

void outputConcentrations(FILE* outFile, double*** electronConcentration, double*** protonConcentration, double*** chargeDensity, int xnumber, int ynumber, int znumber, double plasma_period, double gyroradius) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fprintf(outFile, "%g %g %g\n", electronConcentration[i][j][k]/cube(gyroradius), protonConcentration[i][j][k]/cube(gyroradius), chargeDensity[i][j][k]/(sqrt(cube(gyroradius))*plasma_period));
			}
		}
	}
}

void outputArrayParameter(FILE* outFile, double*** arrayParameter, int xnumber, int ynumber, int znumber) {
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				fprintf(outFile, "%g\n", arrayParameter[i][j][k]);
			}
		}
	}
}

void outputGeneral(FILE* outFile, Simulation* simulation) {
	fprintf(outFile, "%g %g %g %g %g %g %g %g %g\n", simulation->time, simulation->time*simulation->plasma_period, simulation->particleEnergy,
		simulation->electricFieldEnergy, simulation->magneticFieldEnergy, simulation->energy, simulation->momentum.x, simulation->momentum.y, simulation->momentum.z);
}

void outputDivergenceError(FILE* outFile, Simulation* simulation) {
	for(int i = 0; i < simulation->xnumber; ++i) {
		for(int j = 0; j < simulation->ynumber; ++j) {
			for(int k = 0; k < simulation->znumber; ++k) {
				//fprintf(outFile, "%g %g %g\n", div, div - 4*pi*simulation->chargeDensity[i][j][k], div2 - 4*pi*simulation->electricDensity[i][j][k]);
			}
		}
	}
}