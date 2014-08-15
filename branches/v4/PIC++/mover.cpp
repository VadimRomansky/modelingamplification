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

	newparticle.coordinates = newparticle.coordinates + (newparticle.momentum*(deltaT/newparticle.mass));
	newparticle.setMomentumByV(newparticle.velocity() + (oldE + newparticle.velocity().vectorMult(oldB))*beta);

	Vector3d oldV = particle->velocity();
	Vector3d tempV = newparticle.velocity();
	double oldCoordinates[6] = {particle->coordinates.x, particle->coordinates.y, particle->coordinates.z, oldV.x, oldV.y, oldV.z};
	double tempCoordinates[6] = {newparticle.coordinates.x, newparticle.coordinates.y, newparticle.coordinates.z, tempV.x, tempV.y, tempV.z};
	double newCoordinates[6];
	for(int i = 0; i < 6; ++i){
		newCoordinates[i] = oldCoordinates[i];
	}

	double rightPart[6];

	int iterationCount = 0;
	double error = particle->dx/1000;
	while( coordinateDifference(tempCoordinates, newCoordinates, deltaT) > error && iterationCount < 100){
		for(int i = 0; i < 6; ++i){
			tempCoordinates[i] = newCoordinates[i];
		}
		moveParticleNewtonIteration(particle, oldCoordinates, tempCoordinates, newCoordinates);
		iterationCount++;
	}

	if(coordinateDifference(tempCoordinates, newCoordinates, deltaT) > error){
		printf("ERROR newton method did not converge\n");
		newparticle.coordinates = newparticle.coordinates + (newparticle.momentum*(deltaT/newparticle.mass));
		newparticle.setMomentumByV(newparticle.velocity() + (oldE + newparticle.velocity().vectorMult(oldB))*beta);
	}

	*particle = newparticle;
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
	Matrix3d rotationTensor = evaluateAlphaRotationTensor(beta, B);

	tempparticle.coordinates.x += dx;

	Vector3d EderX = (correlationTempEfield(tempparticle) - E)/dx;
	Vector3d BderX = (correlationBfield(tempparticle) - B)/dx;
	Matrix3d rotationTensorDerX = (evaluateAlphaRotationTensor(beta, B) - rotationTensor)/dx;

	tempparticle.coordinates.x -= dx;
	tempparticle.coordinates.y += dy;

	Vector3d EderY = (correlationTempEfield(tempparticle) - E)/dy;
	Vector3d BderY = (correlationBfield(tempparticle) - B)/dy;
	Matrix3d rotationTensorDerY = (evaluateAlphaRotationTensor(beta, B) - rotationTensor)/dy;

	tempparticle.coordinates.y -= dy;
	tempparticle.coordinates.z += dz;

	Vector3d EderZ = (correlationTempEfield(tempparticle) - E)/dz;
	Vector3d BderZ = (correlationBfield(tempparticle) - B)/dz;
	Matrix3d rotationTensorDerZ = (evaluateAlphaRotationTensor(beta, B) - rotationTensor)/dz;

	tempparticle.coordinates.z -= dz;


	Vector3d middleVelocity = rotationTensor*(velocity + E*beta);
	Vector3d middleVelocityDerX = rotationTensorDerX*velocity + (rotationTensorDerX*E + rotationTensor*EderX)*beta;
	Vector3d middleVelocityDerY = rotationTensorDerY*velocity + (rotationTensorDerY*E + rotationTensor*EderY)*beta;
	Vector3d middleVelocityDerZ = rotationTensorDerZ*velocity + (rotationTensorDerZ*E + rotationTensor*EderZ)*beta;
	Vector3d acceleration = (E + (middleVelocity.vectorMult(B))/speed_of_light)*beta;

	functionNewtonMethod[0] = tempCoordinates[0] - oldCoordinates[0] - middleVelocity.x*deltaT;
	functionNewtonMethod[1] = tempCoordinates[1] - oldCoordinates[1] - middleVelocity.y*deltaT;
	functionNewtonMethod[2] = tempCoordinates[2] - oldCoordinates[2] - middleVelocity.z*deltaT;
	functionNewtonMethod[3] = tempCoordinates[3] - oldCoordinates[3] - acceleration.x;
	functionNewtonMethod[4] = tempCoordinates[4] - oldCoordinates[4] - acceleration.y;
	functionNewtonMethod[5] = tempCoordinates[5] - oldCoordinates[5] - acceleration.z;

	leftHalf[0][0] = 1 - middleVelocityDerX.x*0.5*deltaT;
	leftHalf[0][1] = - middleVelocityDerY.x*0.5*deltaT;
	leftHalf[0][2] = - middleVelocityDerZ.x*0.5*deltaT;

	leftHalf[1][0] = - middleVelocityDerX.y*0.5*deltaT;
	leftHalf[1][1] = 1 - middleVelocityDerY.y*0.5*deltaT;
	leftHalf[1][2] = - middleVelocityDerZ.y*0.5*deltaT;

	leftHalf[2][0] = - middleVelocityDerX.z*0.5*deltaT;
	leftHalf[2][1] = - middleVelocityDerY.z*0.5*deltaT;
	leftHalf[2][2] = 1 - middleVelocityDerZ.z*0.5*deltaT;

	leftHalf[3][0] = -0.5*beta*(EderX + (middleVelocityDerX.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderX)/speed_of_light)).x;
	leftHalf[3][1] = -0.5*beta*(EderY + (middleVelocityDerY.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderY)/speed_of_light)).x;
	leftHalf[3][2] = -0.5*beta*(EderZ + (middleVelocityDerZ.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderZ)/speed_of_light)).x;

	leftHalf[4][0] = -0.5*beta*(EderX + (middleVelocityDerX.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderX)/speed_of_light)).y;
	leftHalf[4][1] = -0.5*beta*(EderY + (middleVelocityDerY.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderY)/speed_of_light)).y;
	leftHalf[4][2] = -0.5*beta*(EderZ + (middleVelocityDerZ.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderZ)/speed_of_light)).y;

	leftHalf[5][0] = -0.5*beta*(EderX + (middleVelocityDerX.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderX)/speed_of_light)).z;
	leftHalf[5][1] = -0.5*beta*(EderY + (middleVelocityDerY.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderY)/speed_of_light)).z;
	leftHalf[5][2] = -0.5*beta*(EderZ + (middleVelocityDerZ.vectorMult(B)/speed_of_light) + (middleVelocity.vectorMult(BderZ)/speed_of_light)).z;

	for(int i = 0; i < 6; ++i){
		rightPart[i] = - functionNewtonMethod[i];
		for(int j = 0; j < 3; ++j){
			rightPart[i] += leftHalf[i][j]*tempCoordinates[j];
		}
		if(i > 3){
			rightPart[i] += tempCoordinates[i];
		}
	}


	solveSpecialMatrix(leftHalf, rightPart, newCoordinates);

	for(int i = 0; i < 6; ++i){
		delete[] leftHalf[i];
	}
}