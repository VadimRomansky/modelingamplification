#include "stdlib.h"
#include "math.h"

#include "constants.h"
#include "particle.h"

Particle::Particle(double m, double x0, double y0, double z0, double px0, double py0, double pz0, double dx0, double dy0, double dz0){
	mass = m;

	x = x0;
	y = y0;
	z = z0;

	px = px0;
	py = py0;
	pz = pz0;

	dx = dx0;
	dy = dy0;
	dz = dz0;
}

double Particle::shapeFunctionX(const double& xvalue){
	double epsilon = fabs((x - xvalue)/dx);
	if(epsilon > 1.0){
		return 0;
	} else {
		return 1 - epsilon;
	}
}

double Particle::shapeFunctionY(const double& yvalue){
	double epsilon = fabs((y - yvalue)/dy);
	if(epsilon > 1.0){
		return 0;
	} else {
		return 1 - epsilon;
	}
}

double Particle::shapeFunctionZ(const double& zvalue){
	double epsilon = fabs((z - zvalue)/dz);
	if(epsilon > 1.0){
		return 0;
	} else {
		return 1 - epsilon;
	}
}

double Particle::shapeFunction(const double& xvalue, const double& yvalue, const double& zvalue){
	return shapeFunctionX(xvalue)*shapeFunctionY(yvalue)*shapeFunctionZ(zvalue);
}

double Particle::momentum(){
	return sqrt(px*px + py*py + pz*pz);
}

double Particle::velocityX(){
	double p2 = px*px + py*py + pz*pz;
	double mc2 = mass*speed_of_light*speed_of_light;
	double gamma_factor = sqrt(p2*speed_of_light*speed_of_light + mc2*mc2)/mc2;
	return px/(mass*gamma);
}

double Particle::velocityY(){
	double p2 = px*px + py*py + pz*pz;
	double mc2 = mass*speed_of_light*speed_of_light;
	double gamma_factor = sqrt(p2*speed_of_light*speed_of_light + mc2*mc2)/mc2;
	return py/(mass*gamma);
}
double Particle::velocityZ(){
	double p2 = px*px + py*py + pz*pz;
	double mc2 = mass*speed_of_light*speed_of_light;
	double gamma_factor = sqrt(p2*speed_of_light*speed_of_light + mc2*mc2)/mc2;
	return pz/(mass*gamma);
}