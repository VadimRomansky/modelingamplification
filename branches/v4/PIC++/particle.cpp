#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "particle.h"

Particle::Particle(double m, double q, double w, ParticleTypes t, double x0, double y0, double z0, double px0, double py0, double pz0, double dx0, double dy0, double dz0){
	mass = m;
	charge = q;
	weight = w;
	type = t;

	coordinates.x = x0;
	coordinates.y = y0;
	coordinates.z = z0;

	momentum.x = px0;
	momentum.y = py0;
	momentum.z = pz0;

	dx = dx0;
	dy = dy0;
	dz = dz0;
}

Particle::Particle(const Particle& particle){
	mass = particle.mass;
	charge = particle.charge;
	weight = particle.weight;
	type = particle.type;

	coordinates = particle.coordinates;
	momentum = particle.momentum;

	dx = particle.dx;
	dy = particle.dy;
	dz = particle.dz;
}

double Particle::shapeFunctionX(const double& xvalue){
	double epsilon = fabs((coordinates.x - xvalue)/dx);
	if(epsilon > 1.0){
		return 0;
	} else {
		return 1 - epsilon;
	}
}

double Particle::shapeFunctionY(const double& yvalue){
	double epsilon = fabs((coordinates.y - yvalue)/dy);
	if(epsilon > 1.0){
		return 0;
	} else {
		return 1 - epsilon;
	}
}

double Particle::shapeFunctionZ(const double& zvalue){
	double epsilon = fabs((coordinates.z - zvalue)/dz);
	if(epsilon > 1.0){
		return 0;
	} else {
		return 1 - epsilon;
	}
}

double Particle::shapeFunction(const double& xvalue, const double& yvalue, const double& zvalue){
	return shapeFunctionX(xvalue)*shapeFunctionY(yvalue)*shapeFunctionZ(zvalue);
}

double Particle::momentumAbs(){
	return momentum.norm();
}

Vector3d Particle::velocity(double c){
	double p2 = momentum.x*momentum.x + momentum.y*momentum.y + momentum.z*momentum.z;
	double mc2 = mass*c*c;
	double gamma_factor = sqrt(p2*c*c + mc2*mc2)/mc2;

	return momentum/(mass*gamma_factor);
}

double Particle::velocityX(double c){
	double p2 = momentum.x*momentum.x + momentum.y*momentum.y + momentum.z*momentum.z;
	double mc2 = mass*c*c;
	double gamma_factor = sqrt(p2*c*c + mc2*mc2)/mc2;
	return momentum.x/(mass*gamma_factor);
}

double Particle::velocityY(double c){
	double p2 = momentum.x*momentum.x + momentum.y*momentum.y + momentum.z*momentum.z;
	double mc2 = mass*c*c;
	double gamma_factor = sqrt(p2*c*c + mc2*mc2)/mc2;
	return momentum.y/(mass*gamma_factor);
}
double Particle::velocityZ(double c){
	double p2 = momentum.x*momentum.x + momentum.y*momentum.y + momentum.z*momentum.z;
	double mc2 = mass*c*c;
	double gamma_factor = sqrt(p2*c*c + mc2*mc2)/mc2;
	return momentum.z/(mass*gamma_factor);
}

void Particle::setMomentumByV(Vector3d v, double c){
	if(v.norm() > c){
		printf("ERROR v > c\n");
		exit(0);
	}
	double gamma_factor = 1/sqrt(1 - v.scalarMult(v)/c);
	momentum = v*(mass*gamma_factor);
}

double Particle::gammaFactor(double c){
	double p2 = momentum.x*momentum.x + momentum.y*momentum.y + momentum.z*momentum.z;
	double mc2 = mass*c*c;
	return sqrt(p2*c*c + mc2*mc2)/mc2;
}