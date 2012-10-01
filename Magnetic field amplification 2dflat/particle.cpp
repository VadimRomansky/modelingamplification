#include "stdafx.h"
#include "particle.h"
#include "random.h"
#include "constants.h"
//#include "matrix3d.h"
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
	absoluteY = p.absoluteY;

	absoluteMomentum = p.absoluteMomentum;
	absoluteMomentumX = p.absoluteMomentumX;
	absoluteMomentumY = p.absoluteMomentumY;

	localMomentum = p.localMomentum;
	initialLocalMomentum = p.initialLocalMomentum;
	localMomentumX = p.localMomentumX;
	localMomentumY = p.localMomentumY;

	isCosmicRay = p.isCosmicRay;
	initialMomentum = p.initialMomentum;
	weight = p.weight;
	path = std::list<double>(p.path);
	writePath = p.writePath;
}

Particle::Particle(int a, int znumber, SpaceBin* bin, bool wpath, int n){
	number = n;
	
	double x = bin->x1 + (bin->x2 - bin->x1)*uniRandom();
	double y = bin->y1 + (bin->y2 - bin->y1)*uniRandom();

	absoluteX = x;
	absoluteY = y;
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
	double cosLocalPhi = 2*(uniRandom() - 0.5);
	localMomentumY = localMomentum*sqrt(1- cosLocalTheta*cosLocalTheta)*cosLocalPhi;
	if(localMomentumY != localMomentumY){
		printf("aaa");
	}


	isCosmicRay = false;
	setAbsoluteMomentum(bin->Ux, bin->Uy);
	initialMomentum = absoluteMomentum;
	previousAbsoluteMomentum = initialMomentum;
	path = std::list<double>();
	writePath = wpath;

}

void Particle::setAbsoluteMomentum(double Ux, double Uy){

	double U = sqrt(Ux*Ux + Uy*Uy);

	if(U > epsilon){
		double sqrm = mass*mass;
		double sqrc = c2;
		double sqrp = localMomentum*localMomentum;
		double Vlocal = sqrt(sqrp/(sqrm+sqrp/sqrc));

		double VlocalX = Vlocal*localMomentumX/localMomentum;
		double VlocalY = Vlocal*localMomentumY/localMomentum;
		double VlocalZ;
		if(VlocalX*VlocalX + VlocalY*VlocalY > Vlocal*Vlocal){
			VlocalZ = 0;
		} else {
			VlocalZ = sqrt(Vlocal*Vlocal - VlocalX*VlocalX - VlocalY*VlocalY);
		}

		double gamma = 1/sqrt(1 - U*U/sqrc);
		double denominator = 1 + VlocalX*U/sqrc;

		double VabsoluteRotateX = (VlocalX + U)/denominator;
		double VabsoluteRotateY = VlocalY/(gamma*denominator);

		double cosPhi = Ux/U;
		double sinPhi = Uy/U;

		double VabsoluteX = VabsoluteRotateX*cosPhi - VabsoluteRotateY*sinPhi;
		double VabsoluteY = VabsoluteRotateY*cosPhi + VabsoluteRotateX*sinPhi;
		double VabsoluteZ = VlocalZ/(gamma*denominator);

		double absoluteV = sqrt(VabsoluteX*VabsoluteX + VabsoluteY*VabsoluteY + VabsoluteZ*VabsoluteZ);

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

		if(absoluteV < epsilon){
			printf("aaa");
		}

		absoluteMomentumX = absoluteMomentum*VabsoluteX/absoluteV;
		absoluteMomentumY = absoluteMomentum*VabsoluteY/absoluteV;
	} else {
		absoluteMomentumX = localMomentumX;
		absoluteMomentumY = localMomentumY;
		absoluteMomentum = localMomentum;
	}
}

void Particle::setLocalMomentum(double Ux, double Uy){
	double U = sqrt(Ux*Ux + Uy*Uy);

	if(U > epsilon){
		double sqrm = mass*mass;
		double sqrc = c2;
		double sqrp = absoluteMomentum*absoluteMomentum;

		double Vabsolute = sqrt(sqrp/(sqrm+sqrp/sqrc));

		double VabsoluteX = Vabsolute*absoluteMomentumX/absoluteMomentum;
		double VabsoluteY = Vabsolute*absoluteMomentumY/absoluteMomentum;
		double VabsoluteZ;
		if(VabsoluteX*VabsoluteX + VabsoluteY*VabsoluteY > Vabsolute*Vabsolute){
			VabsoluteZ = 0;
		} else {
			VabsoluteZ = sqrt(Vabsolute*Vabsolute - VabsoluteX*VabsoluteX - VabsoluteY*VabsoluteY);
		}

		double cosPhi = Ux/U;
		double sinPhi = Uy/U;

		double VabsoluteRotateX = VabsoluteX*cosPhi + VabsoluteY*sinPhi;
		double VabsoluteRotateY = VabsoluteY*cosPhi - VabsoluteX*sinPhi;

		double gamma = 1/sqrt(1 - U*U/sqrc);
		double denominator = 1 - VabsoluteRotateX*U/sqrc;

		double VlocalX = (VabsoluteRotateX - U)/denominator;
		double VlocalY = VabsoluteRotateY/(gamma*denominator);

		double VlocalZ = VabsoluteZ/(gamma*denominator);

		double localV = sqrt(VlocalX*VlocalX + VlocalY*VlocalY + VlocalZ*VlocalZ);

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
		if(localV < epsilon){
			printf("aaa");
		}
		localMomentumX = localMomentum*VlocalX/localV;
		localMomentumY = localMomentum*VlocalY/localV;
	} else {
		localMomentumX = absoluteMomentumX;
		localMomentumY = absoluteMomentumY;
		localMomentum = absoluteMomentum;
	}
}

void Particle::setAbsoluteMomentum(SpaceBin* bin){
	setAbsoluteMomentum(bin->Ux, bin->Uy);
}

void Particle::setLocalMomentum(SpaceBin* bin){
	setLocalMomentum(bin->Ux, bin->Uy);
}

double Particle::getAbsoluteV(){
	double sqrm = mass*mass;
	double sqrc = c2;
	double sqrp = absoluteMomentum*absoluteMomentum;
	double v = sqrt(sqrp/(sqrm+sqrp/sqrc));
	if( v != v){
		printf("aaa");
	}
	return v;
}

double Particle::getAbsoluteVX(){
	double sqrm = mass*mass;
	double sqrc = c2;
	double sqrp = absoluteMomentum*absoluteMomentum;
	double v = sqrt(sqrp/(sqrm+sqrp/sqrc));
	if( v != v){
		printf("aaa");
	}
	return v*absoluteMomentumX/absoluteMomentum;
}

double Particle::getAbsoluteVY(){
	double sqrm = mass*mass;
	double sqrc = c2;
	double sqrp = absoluteMomentum*absoluteMomentum;
	double v = sqrt(sqrp/(sqrm+sqrp/sqrc));
	if( v != v){
		printf("aaa");
	}
	return v*absoluteMomentumY/absoluteMomentum;
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
	double sqrc = c2;
	double sqrp = localMomentum*localMomentum;
	double v = sqrt(sqrp/(sqrm+sqrp/sqrc));
	if( v != v){
		printf("v!=v\n");
	}
	return v*localMomentumX/localMomentum;
}

double Particle::getLocalVY(){
	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = localMomentum*localMomentum;
	double v = sqrt(sqrp/(sqrm+sqrp/sqrc));
	if( v != v){
		printf("v!=v\n");
	}
	return v*localMomentumY/localMomentum;
}


double Particle::getEnergy(){
	double v = getAbsoluteV();
	//double c2 = speed_of_light*speed_of_light;
	return (mass*c2*(1/(sqrt(1 - (v*v)/c2)) - 1));
}

double Particle::getInitialEnergy(){
	return sqrt(mass*mass*c2*c2 + initialMomentum*initialMomentum*c2) - mass*c2;
}

void Particle::moveToBinLeftX(SpaceBin* bin){
	absoluteX = bin->x1 + (bin->x2 - bin->x1)*(epsilon);
}

void Particle::moveToBinLeftY(SpaceBin* bin){
	absoluteY = bin->y1 + (bin->y2 - bin->y1)*(epsilon);
}