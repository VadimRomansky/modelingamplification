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
	//for(int i = 0; i < 2; ++i){
		if(i % 100000 == 0) {
			printf("particle number %d\n", i);
		}
		moveParticle(particles[i]);
	}
}

void Simulation::moveParticle(Particle* particle){
	Vector3d E = correlationEfield(particle);
	Vector3d B = correlationBfield(particle);

	double gamma = particle->gammaFactor(speed_of_light_normalized);
	Vector3d velocity = particle->velocity(speed_of_light_normalized);

	particle->oldCoordinates = particle->coordinates;
	particle->coordinates += velocity*deltaT;
	Matrix3d matrix;

	double beta = particle->charge*deltaT/particle->mass;

	matrix.matrix[0][0] = 1;
	matrix.matrix[0][1] = - beta*B.z/(2*gamma*speed_of_light_normalized);
	matrix.matrix[0][2] = beta*B.y/(2*gamma*speed_of_light_normalized);

	matrix.matrix[1][0] = beta*B.z/(2*gamma*speed_of_light_normalized);
	matrix.matrix[1][1] = 1;
	matrix.matrix[1][2] = - beta*B.x/(2*gamma*speed_of_light_normalized);

	matrix.matrix[2][0] = - beta*B.y/(2*gamma*speed_of_light_normalized);
	matrix.matrix[2][1] = beta*B.x/(2*gamma*speed_of_light_normalized);
	matrix.matrix[2][2] = 1;

	Matrix3d inverseMatrix = matrix.inverse();
	Vector3d deltaMomentum = (E + (particle->momentum.vectorMult(B)/(2*gamma*speed_of_light_normalized)))*beta;

	Vector3d rightPart = particle->momentum + (E + (particle->momentum.vectorMult(B)/(2*particle->mass*gamma*speed_of_light_normalized)))*deltaT*particle->charge;

	particle->momentum = matrix.inverse()*rightPart;

	correctParticlePosition(particle);
}

void Simulation::correctParticlePosition(Particle* particle) {
	if(boundaryConditionType == BoundaryConditionTypes::SUPERCONDUCTERLEFT){
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
	} else if(boundaryConditionType == BoundaryConditionTypes::PERIODIC) {
		if(particle->coordinates.x < 0) {
			particle->coordinates.x += xsize;
		}

		if(particle->coordinates.x > xsize) {
			particle->coordinates.x -= xsize;
		}
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