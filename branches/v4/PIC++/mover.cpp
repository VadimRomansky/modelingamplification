#include "stdlib.h"
#include "stdio.h"

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"

void Simulation::moveParticles(){
	printf("moving particles\n");
	//for(int i = 0; i < particles.size(); ++i){
	for(int i = 0; i < 1; ++i){
		moveParticle(particles[i]);
	}
}

void Simulation::moveParticle(Particle* particle){
	Vector3d oldE = correlationTempEfield(particle);
	Vector3d oldB = correlationBfield(particle);

	Particle newparticle = Particle(*particle);

	double beta = particle->charge*deltaT/particle->mass;
	Vector3d oldV = particle->velocity();

	newparticle.coordinates = newparticle.coordinates + oldV*deltaT;
	Vector3d newV = oldV + (oldE + (oldV.vectorMult(oldB))/speed_of_light)*beta;
	newparticle.momentum = newparticle.momentum + (oldE + (oldV.vectorMult(oldB))/speed_of_light)*newparticle.charge*deltaT;

	Vector3d tempV = newparticle.velocity();
	double oldCoordinates[6] = {particle->coordinates.x, particle->coordinates.y, particle->coordinates.z, particle->momentum.x, particle->momentum.y, particle->momentum.z};
	double tempCoordinates[6] = {newparticle.coordinates.x, newparticle.coordinates.y, newparticle.coordinates.z, newparticle.momentum.x, newparticle.momentum.y, newparticle.momentum.z};
	double newCoordinates[6];
	for(int i = 0; i < 6; ++i){
		newCoordinates[i] = tempCoordinates[i];
	}
	moveParticleNewtonIteration(particle, oldCoordinates, tempCoordinates, newCoordinates);

	int iterationCount = 0;
	double error = particle->dx*maxErrorLevel;
	while( coordinateDifference(tempCoordinates, newCoordinates, deltaT, particle->mass) > error && iterationCount < maxNewtonIterations){
		for(int i = 0; i < 6; ++i){
			tempCoordinates[i] = newCoordinates[i];
		}
		moveParticleNewtonIteration(particle, oldCoordinates, tempCoordinates, newCoordinates);
		iterationCount++;
	}

	if(coordinateDifference(tempCoordinates, newCoordinates, deltaT, particle->mass) > error){
	//if(true){
		printf("ERROR newton method did not converge\n");
		//*particle = newparticle;
		//return;

	}

	particle->coordinates.x = newCoordinates[0];
	particle->coordinates.y = newCoordinates[1];
	particle->coordinates.z = newCoordinates[2];
	particle->momentum.x = newCoordinates[3];
	particle->momentum.y = newCoordinates[4];
	particle->momentum.z = newCoordinates[5];

}

void Simulation::moveParticleNewtonIteration(Particle* particle, double* const oldCoordinates, double* const tempCoordinates, double* const newCoordinates){
	double* leftHalf[6];
	double rightPart[6];
	double functionNewtonMethod[6];
	Particle tempparticle = *particle;
	double beta = particle->charge*deltaT/particle->mass;
	double dx = particle->dx/1000;
	double dy = particle->dy/1000;
	double dz = particle->dz/1000;

	for(int i = 0; i < 6; ++i){
		leftHalf[i] = new double[3];
	}

	Vector3d velocity = particle->velocity();

	tempparticle.coordinates.x = (oldCoordinates[0] + tempCoordinates[0])/2;
	tempparticle.coordinates.y = (oldCoordinates[1] + tempCoordinates[1])/2;
	tempparticle.coordinates.z = (oldCoordinates[2] + tempCoordinates[2])/2;

	Vector3d E = correlationTempEfield(tempparticle);
	Vector3d B = correlationBfield(tempparticle);

	double gamma_factor = 1/sqrt(1 - sqr(velocity.getNorm()/speed_of_light));
	double G = (beta*(E.scalarMult(velocity))/speed_of_light_sqr) + gamma_factor;
	double beta1 = beta/G;

	Matrix3d rotationTensor = evaluateAlphaRotationTensor(beta, velocity, E, B);

	Vector3d tempE;
	Vector3d tempB;

	tempparticle.coordinates.x += dx;

	tempE = correlationTempEfield(tempparticle);
	tempB = correlationBfield(tempparticle);
	Vector3d EderX = (tempE - E)/dx;
	Vector3d BderX = (tempB - B)/dx;
	Matrix3d rotationTensorDerX = (evaluateAlphaRotationTensor(beta, velocity, tempE, tempB) - rotationTensor)/dx;

	tempparticle.coordinates.x -= dx;
	tempparticle.coordinates.y += dy;

	tempE = correlationTempEfield(tempparticle);
	tempB = correlationBfield(tempparticle);
	Vector3d EderY = (tempE - E)/dy;
	Vector3d BderY = (tempB - B)/dy;
	Matrix3d rotationTensorDerY = (evaluateAlphaRotationTensor(beta, velocity, tempE, tempB) - rotationTensor)/dy;

	tempparticle.coordinates.y -= dy;
	tempparticle.coordinates.z += dz;

	tempE = correlationTempEfield(tempparticle);
	tempB = correlationBfield(tempparticle);
	Vector3d EderZ = (tempE - E)/dz;
	Vector3d BderZ = (tempB - B)/dz;
	Matrix3d rotationTensorDerZ = (evaluateAlphaRotationTensor(beta, velocity, tempE, tempB) - rotationTensor)/dz;

	tempparticle.coordinates.z -= dz;


	Vector3d middleVelocity = rotationTensor*(velocity*gamma_factor + E*beta1);
	Vector3d middleVelocityDerX = rotationTensorDerX*velocity*gamma_factor + (rotationTensorDerX*E + rotationTensor*EderX)*beta1;
	Vector3d middleVelocityDerY = rotationTensorDerY*velocity*gamma_factor + (rotationTensorDerY*E + rotationTensor*EderY)*beta1;
	Vector3d middleVelocityDerZ = rotationTensorDerZ*velocity*gamma_factor + (rotationTensorDerZ*E + rotationTensor*EderZ)*beta1;
	Vector3d acceleration = (E + (middleVelocity.vectorMult(B))/speed_of_light)*beta;

	functionNewtonMethod[0] = tempCoordinates[0] - oldCoordinates[0] - middleVelocity.x*deltaT;
	functionNewtonMethod[1] = tempCoordinates[1] - oldCoordinates[1] - middleVelocity.y*deltaT;
	functionNewtonMethod[2] = tempCoordinates[2] - oldCoordinates[2] - middleVelocity.z*deltaT;
	functionNewtonMethod[3] = tempCoordinates[3] - oldCoordinates[3] - particle->mass*acceleration.x;
	functionNewtonMethod[4] = tempCoordinates[4] - oldCoordinates[4] - particle->mass*acceleration.y;
	functionNewtonMethod[5] = tempCoordinates[5] - oldCoordinates[5] - particle->mass*acceleration.z;

	leftHalf[0][0] = 1 - middleVelocityDerX.x*0.5*deltaT;
	leftHalf[0][1] = - middleVelocityDerY.x*0.5*deltaT;
	leftHalf[0][2] = - middleVelocityDerZ.x*0.5*deltaT;

	leftHalf[1][0] = - middleVelocityDerX.y*0.5*deltaT;
	leftHalf[1][1] = 1 - middleVelocityDerY.y*0.5*deltaT;
	leftHalf[1][2] = - middleVelocityDerZ.y*0.5*deltaT;

	leftHalf[2][0] = - middleVelocityDerX.z*0.5*deltaT;
	leftHalf[2][1] = - middleVelocityDerY.z*0.5*deltaT;
	leftHalf[2][2] = 1 - middleVelocityDerZ.z*0.5*deltaT;

	leftHalf[3][0] = -0.5*particle->mass*beta*(EderX + (middleVelocityDerX.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderX)/speed_of_light)).x;
	leftHalf[3][1] = -0.5*particle->mass*beta*(EderY + (middleVelocityDerY.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderY)/speed_of_light)).x;
	leftHalf[3][2] = -0.5*particle->mass*beta*(EderZ + (middleVelocityDerZ.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderZ)/speed_of_light)).x;

	leftHalf[4][0] = -0.5*particle->mass*beta*(EderX + (middleVelocityDerX.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderX)/speed_of_light)).y;
	leftHalf[4][1] = -0.5*particle->mass*beta*(EderY + (middleVelocityDerY.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderY)/speed_of_light)).y;
	leftHalf[4][2] = -0.5*particle->mass*beta*(EderZ + (middleVelocityDerZ.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderZ)/speed_of_light)).y;

	leftHalf[5][0] = -0.5*particle->mass*beta*(EderX + (middleVelocityDerX.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderX)/speed_of_light)).z;
	leftHalf[5][1] = -0.5*particle->mass*beta*(EderY + (middleVelocityDerY.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderY)/speed_of_light)).z;
	leftHalf[5][2] = -0.5*particle->mass*beta*(EderZ + (middleVelocityDerZ.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderZ)/speed_of_light)).z;

	for(int i = 0; i < 6; ++i){
		rightPart[i] = - functionNewtonMethod[i];
		for(int j = 0; j < 3; ++j){
			rightPart[i] += leftHalf[i][j]*tempCoordinates[j];
		}
		if(i >= 3){
			rightPart[i] += tempCoordinates[i];
		}
	}


	solveSpecialMatrix(leftHalf, rightPart, newCoordinates);

	for(int i = 0; i < 6; ++i){
		delete[] leftHalf[i];
	}
}