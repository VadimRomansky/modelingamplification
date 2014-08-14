#include "stdlib.h"
#include "stdio.h"

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"

void Simulation::moveParticles(){
	for(int i = 0; i < particles.size(); ++i){
		moveParticle(particles[i]);
	}
}

void Simulation::moveParticle(Particle* particle){
	Vector3d E = correlationTempEfield(particle);
	Vector3d B = correlationBfield(particle);

	Particle* newparticle = new Particle(*particle);

	double beta = particle->charge*deltaT/particle->mass;

	Matrix3d rotationTensor = evaluateAlphaRotationTensor(beta,B);

	newparticle->coordinates = newparticle->coordinates + (newparticle->momentum*(deltaT/newparticle->mass));
	newparticle->setMomentumByV(newparticle->velocity() + (E + newparticle->velocity().vectorMult(B))*beta);
}