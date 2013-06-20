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
	absoluteY = p.absoluteY;
	absoluteZ = p.absoluteZ;

	absoluteMomentum = p.absoluteMomentum;
	absoluteMomentumTheta = p.absoluteMomentumTheta;
	absoluteMomentumPhi = p.absoluteMomentumPhi;

	localMomentum = p.localMomentum;
	initialLocalMomentum = p.initialLocalMomentum;
	localMomentumX = p.localMomentumX;
	localMomentumY = p.localMomentumY;
	localMomentumZ = p.localMomentumZ;

	isCosmicRay = p.isCosmicRay;
	initialMomentum = p.initialMomentum;
	weight = p.weight;
	path = std::list<double>(p.path);
	writePath = p.writePath;
}

Particle::Particle(double r, double temperature, int a, int znumber){
	double x = uniRandom() - 0.5;
	double y = uniRandom() - 0.5;
	double z = uniRandom() - 0.5;
	
	double theta = acos(z/sqrt(x*x + y*y + z*z));
	if(abs(x*x + y*y + z*z) < epsilon){
		theta = pi/2;
	}
	double phi = atan2(y,x);

	absoluteX = r*sin(theta)*cos(phi);
	absoluteY = r*sin(theta)*sin(phi);
	absoluteZ = r*cos(theta);

	A = a;
	Z = znumber;
	mass = A*massProton;
	weight = 1;
	double px = randomMaxwell(temperature, A*(massProton));
	double py = randomMaxwell(temperature, A*(massProton));
	double pz = randomMaxwell(temperature, A*(massProton));

	localMomentum = sqrt(px*px+py*py+pz*pz);
	initialLocalMomentum = localMomentum;
	localMomentumX = px;
	localMomentumY = py;
	localMomentumZ = pz;

	isCosmicRay = false;
	path = std::list<double>();
	writePath = false;
}


Particle::Particle(double x, double y, double z, double temperature, int a, int znumber){
	absoluteX = x;
	absoluteY = y;
	absoluteZ = z;

	A = a;
	Z = znumber;
	mass = A*massProton;
	weight = 1;
	double px = randomMaxwell(temperature, A*(massProton));
	double py = randomMaxwell(temperature, A*(massProton));
	double pz = randomMaxwell(temperature, A*(massProton));

	localMomentum = sqrt(px*px+py*py+pz*pz);
	initialLocalMomentum = localMomentum;
	localMomentumX = px;
	localMomentumY = py;
	localMomentumZ = pz;

	isCosmicRay = false;
	path = std::list<double>();
	writePath = false;
}


Particle::Particle(double x, double y, double z, double temperature, int a, int znumber, double Ur, double Utheta, double Uphi){
	absoluteX = x;
	absoluteY = y;
	absoluteZ = z;
	A = a;
	Z = znumber;
	mass = A*massProton;
	weight = 1;
	double px = randomMaxwell(temperature, A*(massProton));
	double py = randomMaxwell(temperature, A*(massProton));
	double pz = randomMaxwell(temperature, A*(massProton));

	localMomentum = sqrt(px*px+py*py+pz*pz);
	initialLocalMomentum = localMomentum;
	localMomentumX = px;
	localMomentumY = py;
	localMomentumZ = pz;
	isCosmicRay = false;
	setAbsoluteMomentum(Ur,Utheta,Uphi);
	initialMomentum = absoluteMomentum;
	previousAbsoluteMomentum = initialMomentum;
	path = std::list<double>();
	writePath = false;
}

Particle::Particle(double x, double y, double z, double temperature, int a, int znumber, double Ur, double Utheta, double Uphi, bool wPath){
	absoluteX = x;
	absoluteY = y;
	absoluteZ = z;
	A = a;
	Z = znumber;
	mass = A*massProton;
	weight = 1;
	double px = randomMaxwell(temperature, A*(massProton));
	double py = randomMaxwell(temperature, A*(massProton));
	double pz = randomMaxwell(temperature, A*(massProton));

	localMomentum = sqrt(px*px+py*py+pz*pz);
	initialLocalMomentum = localMomentum;
	localMomentumX = px;
	localMomentumY = py;
	localMomentumZ = pz;
	isCosmicRay = false;
	setAbsoluteMomentum(Ur,Utheta,Uphi);
	initialMomentum = absoluteMomentum;
	path = std::list<double>();
	writePath = wPath;
}

Particle::Particle(double x, double y, double z, double temperature, int a, int znumber, double Ur, double Utheta, double Uphi, double px,double py, double pz, double ptheta, double pphi){
	absoluteX = x;
	absoluteY = y;
	absoluteZ = z;
	A = a;
	Z = znumber;
	mass = A*massProton;
	weight = 1;
	localMomentum = sqrt(px*px + py*py + pz*pz);
	initialLocalMomentum = localMomentum;
	localMomentumX = px;
	localMomentumY = py;
	localMomentumZ = pz;

	isCosmicRay = false;
	setAbsoluteMomentum(Ur,Utheta,Uphi);
	initialMomentum = absoluteMomentum;
	previousAbsoluteMomentum = initialMomentum;
	path = std::list<double>();
	writePath = false;
}

Particle::Particle(double x, double y, double z, double temperature, int a, int znumber, double U, double Utheta, double Uphi, double px, double py, double pz, bool wPath){
	absoluteX = x;
	absoluteY = y;
	absoluteZ = z;
	A = a;
	Z = znumber;
	mass = A*massProton;
	weight = 1;
	localMomentum = sqrt(px*px + py*py + pz*pz);
	initialLocalMomentum = localMomentum;
	localMomentumX = px;
	localMomentumY = py;
	localMomentumZ = pz;

	isCosmicRay = false;
	setAbsoluteMomentum(U,Utheta,Uphi);
	initialMomentum = absoluteMomentum;
	previousAbsoluteMomentum = initialMomentum;
	path = std::list<double>();
	writePath = wPath;
}

void Particle::setAbsoluteMomentum(double U, double Utheta, double Uphi){
	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = localMomentum*localMomentum;
	double V = sqrt(sqrp/(sqrm+sqrp/sqrc));
	//double c2 = speed_of_light*speed_of_light;
	double theta = acos(localMomentumZ/localMomentum);
	if(abs(localMomentum) < epsilon*sqrt((kBoltzman)*1000*massProton)){
		theta = pi/2;
	}
	double phi = atan2(localMomentumY,localMomentumX);

	double ux;
	double uy;
	double uz;

	if((thetagridNumber == 1)&&(phigridNumber == 1)){
		ux = U*sin(getAbsoluteTheta())*cos(getAbsolutePhi());
		uy = U*sin(getAbsoluteTheta())*sin(getAbsolutePhi());
		uz = U*cos(getAbsoluteTheta());		
	} else {
		//double c2 = speed_of_light*speed_of_light;

		ux = U*sin(Utheta)*cos(Uphi);
		uy = U*sin(Utheta)*sin(Uphi);
		uz = U*cos(Utheta);
	}	
	
	//u - скорость плазмы в абсолютной СО. Локальная ось z направлена по скорости плазмы
	Matrix3d* matrix = Matrix3d::createBasisByOneVector(vector3d(ux,uy,uz));

	vector3d particleV = vector3d(V*sin(theta)*cos(phi),V*sin(theta)*sin(phi),V*cos(theta));

	vector3d particleAbsoluteV = summVelocity(particleV,-U);

	//Matrix3d* invertMatrix = matrix->Inverse();

	double absoluteV = particleAbsoluteV.getNorm();

	particleAbsoluteV = matrix->multiply(particleAbsoluteV);

	absoluteV = particleAbsoluteV.getNorm();

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

	absoluteMomentumTheta = acos(particleAbsoluteV.z/absoluteV);
	if(abs(absoluteV) < epsilon){
		absoluteMomentumTheta = pi/2;
	}
	absoluteMomentumPhi = atan2(particleAbsoluteV.y, particleAbsoluteV.x);
	
	delete matrix;
	//delete invertMatrix;
}

void Particle::setLocalMomentum(double U, double Utheta,double Uphi){
	double sqrm = mass*mass;
	double sqrp = absoluteMomentum*absoluteMomentum;
	double V = sqrt(sqrp/(sqrm+sqrp/c2));
	double theta = absoluteMomentumTheta;
	double phi = absoluteMomentumPhi;
	double ux;
	double uy;
	double uz;
	if((thetagridNumber == 1)&&(phigridNumber == 1)){
		ux = U*sin(getAbsoluteTheta())*cos(getAbsolutePhi());
		uy = U*sin(getAbsoluteTheta())*sin(getAbsolutePhi());
		uz = U*cos(getAbsoluteTheta());		
	} else {
		//double c2 = speed_of_light*speed_of_light;

		ux = U*sin(Utheta)*cos(Uphi);
		uy = U*sin(Utheta)*sin(Uphi);
		uz = U*cos(Utheta);
	}

	double localV1 = getLocalV();

	Matrix3d* matrix = Matrix3d::createBasisByOneVector(vector3d(ux,uy,uz));

	vector3d particleAbsoluteV = vector3d(V*sin(theta)*cos(phi),V*sin(theta)*sin(phi),V*cos(theta));

	Matrix3d* invertMatrix = matrix->Inverse();

	double absoluteV1 = particleAbsoluteV.getNorm();

	particleAbsoluteV = invertMatrix->multiply(particleAbsoluteV);

	double absoluteV = particleAbsoluteV.getNorm();

	vector3d particleLocalV = summVelocity(particleAbsoluteV,U);

	//for debug with U = const only!!!!!!!!!!!!!!
	//if(abs(localMomentum - initialLocalMomentum) > epsilon*localMomentum){
		//printf("aaa");
	//}


	double localV = particleLocalV.getNorm();
	if(localV > speed_of_light){
		if(localV < (1 + epsilon)*speed_of_light){
			localV = (1 - epsilon)*speed_of_light;
		} else {
			printf("localV > c\n");
		}
	}

	localMomentum = mass*localV/sqrt(1 - localV*localV/c2);

	//if(abs(localMomentum - initialLocalMomentum) > epsilon*localMomentum){
		//printf("aaa");
	//}

	if(localMomentum != localMomentum){
		printf("localMomentum != localMomentum\n");
	}
	localMomentumX = mass*particleLocalV.x/sqrt(1 - localV*localV/c2);
	localMomentumY = mass*particleLocalV.y/sqrt(1 - localV*localV/c2);
	localMomentumZ = mass*particleLocalV.z/sqrt(1 - localV*localV/c2);

	delete matrix;
	delete invertMatrix;
}

void Particle::setAbsoluteMomentum(SpaceBin* bin){
	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = localMomentum*localMomentum;
	double V = sqrt(sqrp/(sqrm+sqrp/sqrc));
	//double c2 = speed_of_light*speed_of_light;
	double theta = acos(localMomentumZ/localMomentum);
	if(abs(localMomentum) < epsilon*sqrt((kBoltzman)*1000*massProton)){
		theta = pi/2;
	}
	double phi = atan2(localMomentumY,localMomentumX);

	double ux;
	double uy;
	double uz;

	if((thetagridNumber == 1)&&(phigridNumber == 1)){
		ux = 0;
		uy = 0;
		uz = bin->U;		
	} else {
		//double c2 = speed_of_light*speed_of_light;
	}	
	
	//u - скорость плазмы в абсолютной СО. Локальная ось z направлена по скорости плазмы
	Matrix3d* matrix = bin->matrix;

	vector3d particleV = vector3d(V*sin(theta)*cos(phi),V*sin(theta)*sin(phi),V*cos(theta));

	vector3d particleAbsoluteV = summVelocity(particleV,-bin->U);

	//Matrix3d* invertMatrix = matrix->Inverse();

	double absoluteV = particleAbsoluteV.getNorm();

	particleAbsoluteV = matrix->multiply(particleAbsoluteV);

	absoluteV = particleAbsoluteV.getNorm();

	if(absoluteV > speed_of_light){
		if(absoluteV < (1 + epsilon)*speed_of_light){
			absoluteV = (1 - epsilon)*speed_of_light;
			printf("v = c\n");
		} else {
			printf("v > c\n");
		}
	}

	if(abs(1 - absoluteV*absoluteV/c2) < epsilon){
		printf("absoluteV = c in setAbsoluteMomentum\n");
	}
	absoluteMomentum = mass*absoluteV/sqrt(1 - absoluteV*absoluteV/c2);
	if(absoluteMomentum != absoluteMomentum){
		printf("absoluteMomentum != absoluteMomentum\n");
	}

	absoluteMomentumTheta = acos(particleAbsoluteV.z/absoluteV);
	if(abs(absoluteV) < epsilon){
		absoluteMomentumTheta = pi/2;
	}
	absoluteMomentumPhi = atan2(particleAbsoluteV.y, particleAbsoluteV.x);
	
	//delete matrix;
	//delete invertMatrix;
}

void Particle::setLocalMomentum(SpaceBin* bin){
	double sqrm = mass*mass;
	double sqrp = absoluteMomentum*absoluteMomentum;
	double V = sqrt(sqrp/(sqrm+sqrp/c2));
	double theta = absoluteMomentumTheta;
	double phi = absoluteMomentumPhi;
	double ux;
	double uy;
	double uz;
	if((thetagridNumber == 1)&&(phigridNumber == 1)){
		ux = 0;
		uy = 0;
		uz = bin->U;		
	} else {
		//double c2 = speed_of_light*speed_of_light;
	}

	double localV1 = getLocalV();

	Matrix3d* matrix =bin->matrix;

	vector3d particleAbsoluteV = vector3d(V*sin(theta)*cos(phi),V*sin(theta)*sin(phi),V*cos(theta));

	Matrix3d* invertMatrix = bin->invertMatrix;

	double absoluteV1 = particleAbsoluteV.getNorm();

	particleAbsoluteV = invertMatrix->multiply(particleAbsoluteV);

	double absoluteV = particleAbsoluteV.getNorm();

	vector3d particleLocalV = summVelocity(particleAbsoluteV,bin->U);

	//for debug with U = const only!!!!!!!!!!!!!!
	//if(abs(localMomentum - initialLocalMomentum) > epsilon*localMomentum){
		//printf("aaa");
	//}


	double localV = particleLocalV.getNorm();
	if(localV > speed_of_light){
		if(localV < (1 + epsilon)*speed_of_light){
			localV = (1 - epsilon)*speed_of_light;
		} else {
			printf("localV > c\n");
		}
	}


	if(abs(1 - localV*localV/c2) < epsilon){
		printf("localV = c in setLocalMomentum\n");
	}
	localMomentum = mass*localV/sqrt(1 - localV*localV/c2);

	//if(abs(localMomentum - initialLocalMomentum) > epsilon*localMomentum){
		//printf("aaa");
	//}

	if(localMomentum != localMomentum){
		printf("localMomentum != localMomentum\n");
	}
	localMomentumX = mass*particleLocalV.x/sqrt(1 - localV*localV/c2);
	localMomentumY = mass*particleLocalV.y/sqrt(1 - localV*localV/c2);
	localMomentumZ = mass*particleLocalV.z/sqrt(1 - localV*localV/c2);

	//delete matrix;
	//delete invertMatrix;
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

double Particle::getAbsoluteR(){
	return sqrt(absoluteX*absoluteX + absoluteY*absoluteY + absoluteZ*absoluteZ);
}

double Particle::getAbsoluteTheta(){
	double r = getAbsoluteR();
	if(abs(r) < epsilon){
		return pi/2;
	}
	return acos(absoluteZ/r);
}

double Particle::getAbsolutePhi(){
	double phi = atan2(absoluteY, absoluteX);
	if (phi < 0) {
		phi = phi + 2*pi;
	}
	return phi;
}

double Particle::getEnergy(){
	double v = getAbsoluteV();
	//double c2 = speed_of_light*speed_of_light;
	return (mass*c2/(sqrt(1 - (v*v)/c2)) - mass*c2);
}

double Particle::getRadialSpeed(){
	return getAbsoluteV()*(cos(absoluteMomentumPhi)*sin(absoluteMomentumTheta)*cos(getAbsolutePhi())*sin(getAbsoluteTheta()) + sin(absoluteMomentumPhi)*sin(absoluteMomentumTheta)*sin(getAbsolutePhi())*sin(getAbsoluteTheta()) + cos(absoluteMomentumTheta)*cos(getAbsoluteTheta()));
}

double Particle::getAbsoluteVR(){
	return getAbsoluteV()*(cos(absoluteMomentumPhi)*sin(absoluteMomentumTheta)*cos(getAbsolutePhi())*sin(getAbsoluteTheta()) + sin(absoluteMomentumPhi)*sin(absoluteMomentumTheta)*sin(getAbsolutePhi())*sin(getAbsoluteTheta()) + cos(absoluteMomentumTheta)*cos(getAbsoluteTheta()));
}

double Particle::getAbsoluteVPhi(){
	double absolutePhi = getAbsolutePhi();
	return getAbsoluteV()*(-sin(absolutePhi)*sin(absoluteMomentumTheta)*cos(absoluteMomentumPhi)+cos(absolutePhi)*sin(absoluteMomentumTheta)*sin(absoluteMomentumPhi));
}

double Particle::getAbsoluteVTheta(){
	double absolutePhi = getAbsolutePhi();
	double absoluteTheta = getAbsoluteTheta();
	return getAbsoluteV()*(-sin(absoluteTheta)*cos(absoluteMomentumPhi)+cos(absoluteTheta)*cos(absolutePhi)*sin(absoluteMomentumTheta)*cos(absoluteMomentumPhi)+cos(absoluteTheta)*sin(absolutePhi)*sin(absoluteMomentumTheta)*sin(absoluteMomentumPhi));
}

double Particle::getInitialEnergy(){
	return sqrt(mass*mass*c2*c2 + initialMomentum*initialMomentum*c2) - mass*c2;
}