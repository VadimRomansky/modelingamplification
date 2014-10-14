#include "stdlib.h"
#include "stdio.h"

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"

void Simulation::moveParticles(){
	printf("moving particles\n");
	for(int i = 0; i < particles.size(); ++i){
	//for(int i = 0; i < 1; ++i){
		if(i % 100000 == 0) {
			printf("particle number %d\n", i);
		}
		moveParticle(particles[i]);
	}
}

void Simulation::moveParticle(Particle* particle){
	

	correctParticlePosition(particle);
}

void Simulation::correctParticlePosition(Particle* particle) {
	if(particle->coordinates.x < 0) {
		particle->coordinates.x = -particle->coordinates.x;
		particle->momentum.x = -particle->momentum.x;
	}
	if(particle->coordinates.x > xsize) {
		std::vector<Particle*>::iterator it = particles.begin();

		while(it != particles.end()) {
			if(*it == particle) {
				break;
			}
			++it;
		}
		particles.erase(it);
		particlesNumber--;
		delete particle;
		return;
	}
	if(particle->coordinates.y < 0) {
		particle->coordinates.y += ysize;
	}
	if(particle->coordinates.y > ysize) {
		particle->coordinates.y -= ysize;
	}
	if(particle->coordinates.z < 0) {
		particle->coordinates.z += zsize;
	}
	if(particle->coordinates.z > zsize){
		particle->coordinates.z -= zsize;
	}
}