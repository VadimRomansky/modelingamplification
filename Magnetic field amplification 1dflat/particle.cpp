#include "stdafx.h"
#include "particle.h"
#include "random.h"
#include "constants.h"
#include "matrix3d.h"
#include "SpaceBin.h"
#include "util.h"

Particle::Particle(){
	weight = 1;
}

Particle::Particle(const Particle& p){
	A = p.A;
	Z = p.Z;
	mass = p.mass;
	absoluteX = p.absoluteX;

	absoluteMomentum = p.absoluteMomentum;
	absoluteMomentumX = p.absoluteMomentumX;

	localMomentum = p.localMomentum;
	initialLocalMomentum = p.initialLocalMomentum;
	localMomentumX = p.localMomentumX;

	isCosmicRay = p.isCosmicRay;
	initialMomentum = p.initialMomentum;
	weight = p.weight;
	path = std::list<double>(p.path);
	writePath = p.writePath;
}

Particle::Particle(double x, double temperature, int a, int znumber){

	
	absoluteX = x;

	A = a;
	Z = znumber;
	mass = A*massProton;
	weight = 1;
	double px = randomMaxwell(temperature, A*(massProton));
	double py = randomMaxwell(temperature, A*(massProton));
	double pz = randomMaxwell(temperature, A*(massProton));
	//todo
	double cosLocalTheta = 2*(uniRandom() - 0.5);
	localMomentum = sqrt(px*px+py*py+pz*pz);
	initialLocalMomentum = localMomentum;
	localMomentumX = localMomentum*cosLocalTheta;


	isCosmicRay = false;
	path = std::list<double>();
	writePath = false;
}


Particle::Particle(double x, double temperature, int a, int znumber, double U){
	absoluteX = x;

	A = a;
	Z = znumber;
	mass = A*massProton;
	weight = 1;
	double px = randomMaxwell(temperature, A*(massProton));
	double py = randomMaxwell(temperature, A*(massProton));
	double pz = randomMaxwell(temperature, A*(massProton));

	double cosLocalTheta = 2*(uniRandom() - 0.5);
	localMomentum = sqrt(px*px+py*py+pz*pz);
	initialLocalMomentum = localMomentum;
	localMomentumX = localMomentum*cosLocalTheta;

	isCosmicRay = false;
	setAbsoluteMomentum(U);
	initialMomentum = absoluteMomentum;
	previousAbsoluteMomentum = initialMomentum;
	path = std::list<double>();
	writePath = false;
}

Particle::Particle(double x, double temperature, int a, int znumber, double U, bool wPath){
	absoluteX = x;
	A = a;
	Z = znumber;
	mass = A*massProton;
	weight = 1;
	double px = randomMaxwell(temperature, A*(massProton));
	double py = randomMaxwell(temperature, A*(massProton));
	double pz = randomMaxwell(temperature, A*(massProton));

	double cosLocalTheta = 2*(uniRandom() - 0.5);
	localMomentum = sqrt(px*px+py*py+pz*pz);
	initialLocalMomentum = localMomentum;
	localMomentumX = localMomentum*cosLocalTheta;

	isCosmicRay = false;
	setAbsoluteMomentum(U);
	initialMomentum = absoluteMomentum;
	path = std::list<double>();
	writePath = wPath;
}

Particle::Particle(double x, double temperature, int a, int znumber, double U, double px,double py, double pz, bool wpath){
	absoluteX = x;

	A = a;
	Z = znumber;
	mass = A*massProton;
	weight = 1;
	double cosLocalTheta = 2*(uniRandom() - 0.5);
	localMomentum = sqrt(px*px+py*py+pz*pz);
	initialLocalMomentum = localMomentum;
	localMomentumX = localMomentum*cosLocalTheta;


	isCosmicRay = false;
	setAbsoluteMomentum(U);
	initialMomentum = absoluteMomentum;
	previousAbsoluteMomentum = initialMomentum;
	path = std::list<double>();
	writePath = false;
}

Particle::Particle(int a, int znumber, SpaceBin* bin, bool wpath, int n){
	number = n;
	
	double x = bin->r1 + (bin->r2 - bin->r1)*uniRandom();

	absoluteX = x;
	A = a;
	Z = znumber;
	mass = A*massProton;
	weight = bin->density*bin->volume/mass;

	double px = randomMaxwell(bin->temperature, mass);
	double py = randomMaxwell(bin->temperature, mass);
	double pz = randomMaxwell(bin->temperature, mass);

	double cosLocalTheta = 2*(uniRandom() - 0.5);
	localMomentum = sqrt(px*px+py*py+pz*pz);
	initialLocalMomentum = localMomentum;
	localMomentumX = localMomentum*cosLocalTheta;


	isCosmicRay = false;
	setAbsoluteMomentum(bin->U);
	initialMomentum = absoluteMomentum;
	previousAbsoluteMomentum = initialMomentum;
	path = std::list<double>();
	writePath = wpath;

}

void Particle::setAbsoluteMomentum(double U){
	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = localMomentum*localMomentum;
	double V = sqrt(sqrp/(sqrm+sqrp/sqrc));

	double Vx = V*localMomentumX/localMomentum;
	double Vy;
	if(abs(Vx) >= V){
		Vy = 0;
	} else {
		Vy = sqrt(V*V - Vx*Vx);
	}

	double Ux = (Vx + U)/(1 + Vx*U/sqrc);
	double Uy = Vy*sqrt(1 - U*U/sqrc)/(1 + Vx*U/sqrc);

	double absoluteV = sqrt(Ux*Ux + Uy*Uy);


	
	
	//u - скорость плазмы в абсолютной СО. Локальная ось z направлена по скорости плазмы


	if(absoluteV > speed_of_light){
		if(absoluteV < (1 + epsilon)*speed_of_light){
			absoluteV = (1 - epsilon)*speed_of_light;
			printf("v = c\n");
		} else {
			printf("v > c\n");
		}
	}

	absoluteMomentum = mass*absoluteV/sqrt(1 - absoluteV*absoluteV/c2);
	if(absoluteMomentum != absoluteMomentum){
		printf("absoluteMomentum != absoluteMomentum\n");
	}

	absoluteMomentumX = absoluteMomentum*Ux/absoluteV;

	//delete invertMatrix;
}

void Particle::setLocalMomentum(double U){
	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = absoluteMomentum*absoluteMomentum;

	double V = sqrt(sqrp/(sqrm+sqrp/sqrc));

	double Vx = V*absoluteMomentumX/absoluteMomentum;
	double Vy;
	if(abs(Vx) >= V){
		Vy = 0;
	} else {
		Vy = sqrt(V*V - Vx*Vx);
	}

	double Ux = (Vx - U)/(1 - Vx*U/sqrc);
	double Uy = Vy*sqrt(1 - U*U/sqrc)/(1 - Vx*U/sqrc);

	double localV = sqrt(Ux*Ux + Uy*Uy);


	
	
	//u - скорость плазмы в абсолютной СО. Локальная ось z направлена по скорости плазмы


	if(localV > speed_of_light){
		if(localV < (1 + epsilon)*speed_of_light){
			localV = (1 - epsilon)*speed_of_light;
			printf("v = c\n");
		} else {
			printf("v > c\n");
		}
	}

	localMomentum = mass*localV/sqrt(1 - localV*localV/c2);
	if(localMomentum != localMomentum){
		printf("localMomentum != localMomentum\n");
	}

	localMomentumX = localMomentum*Ux/localV;
}

void Particle::setAbsoluteMomentum(SpaceBin* bin){
	setAbsoluteMomentum(bin->U);
}

void Particle::setLocalMomentum(SpaceBin* bin){
	setLocalMomentum(bin->U);
}

double Particle::getAbsoluteV(){
	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = absoluteMomentum*absoluteMomentum;
	double v = sqrt(sqrp/(sqrm+sqrp/sqrc));
	if( v != v){
		printf("aaa");
	}
	return v;
}

double Particle::getAbsoluteVX(){
	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = absoluteMomentum*absoluteMomentum;
	double v = sqrt(sqrp/(sqrm+sqrp/sqrc));
	if( v != v){
		printf("aaa");
	}
	return v*absoluteMomentumX/absoluteMomentum;
}

double Particle::getLocalV(){
	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = localMomentum*localMomentum;
	double v = sqrt(sqrp/(sqrm+sqrp/sqrc));
	if( v != v){
		printf("v!=v\n");
	}
	return v;
}

double Particle::getLocalVX(){
	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = localMomentum*localMomentum;
	double v = sqrt(sqrp/(sqrm+sqrp/sqrc));
	if( v != v){
		printf("v!=v\n");
	}
	return v*localMomentumX/localMomentum;
}


double Particle::getEnergy(){
	double v = getAbsoluteV();
	//double c2 = speed_of_light*speed_of_light;
	return (mass*c2*(1/(sqrt(1 - (v*v)/c2)) - 1));
}

double Particle::getInitialEnergy(){
	return sqrt(mass*mass*c2*c2 + initialMomentum*initialMomentum*c2) - mass*c2;
}

void Particle::moveToBinRight(SpaceBin* bin){
	absoluteX = bin->r1 + (bin->r2 - bin->r1)*(1-epsilon);
}

void Particle::moveToBinLeft(SpaceBin* bin){
	absoluteX = bin->r1 + (bin->r2 - bin->r1)*(epsilon);
}