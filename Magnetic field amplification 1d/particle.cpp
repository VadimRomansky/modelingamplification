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
	previousAbsoluteMomentum = p.previousAbsoluteMomentum;
	weight = p.weight;
	//path = std::list<double>(p.path);
	writePath = p.writePath;
}

Particle::Particle(double r, double temperature, int a, int znumber){
	double x = uniRandom() - 0.5;
	double y = uniRandom() - 0.5;
	double z = uniRandom() - 0.5;
	
	double theta = acos(z/sqrt(x*x + y*y + z*z));
	if(abs(x*x + y*y + z*z) < DBL_EPSILON){
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
	previousAbsoluteMomentum = initialMomentum;
	localMomentumX = px;
	localMomentumY = py;
	localMomentumZ = pz;

	isCosmicRay = false;
	//path = std::list<double>();
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
	previousAbsoluteMomentum = initialMomentum;
	localMomentumX = px;
	localMomentumY = py;
	localMomentumZ = pz;

	isCosmicRay = false;
	//path = std::list<double>();
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
	//path = std::list<double>();
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
	previousAbsoluteMomentum = initialMomentum;
	//path = std::list<double>();
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
	//path = std::list<double>();
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
	//path = std::list<double>();
	writePath = wPath;
}

Particle::Particle(int a, int znumber, SpaceBin* bin, bool wpath){
	double x = startX;
	double y = startY;
	double z = bin->r1 + (bin->r2 - bin->r1)*uniRandom();
	double theta = acos(z/sqrt(x*x + y*y +z*z));
	if(abs(x*x + y*y + z*z) < DBL_EPSILON){
		theta = pi/2;
	}
	double phi = atan2(y,x);
	if(phi < 0){
		phi = phi + 2*pi;
	}
	double r = sqrt(x*x + y*y + z*z);
	/*while(!(order(bin->r1,r,bin->r2) && order(bin->theta1, theta, bin->theta2) && order(bin->phi1, phi, bin->phi2))){
		x = bin->r2*2*(uniRandom() - 0.5);
		y = bin->r2*2*(uniRandom() - 0.5);
		z = bin->r2*2*(uniRandom() - 0.5);
		theta = acos(z/sqrt(x*x + y*y +z*z));
		phi = atan2(y,x);
		if(phi < 0){
			phi = phi + 2*pi;
		}
		r = sqrt(x*x + y*y + z*z);
	}*/

	absoluteX = x;
	absoluteY = y;
	absoluteZ = z;
	A = a;
	Z = znumber;
	mass = A*massProton;
	weight = bin->density*bin->volume/mass;
	//weight = 1.0;

	double px = randomMaxwell(bin->temperature, mass);
	double py = randomMaxwell(bin->temperature, mass);
	double pz = randomMaxwell(bin->temperature, mass);

	localMomentum = sqrt(px*px+py*py+pz*pz);
	initialLocalMomentum = localMomentum;

	double localMomentumPhi = 2*pi*uniRandom();
	if(localMomentumPhi >= 2*pi){
		localMomentumPhi -=2*pi;
	}
	double localMomentumCosTheta = 2*(uniRandom() - 0.5);
	if( localMomentumCosTheta > 1){
		localMomentumCosTheta = 1;
	}
	if( localMomentumCosTheta < -1){
		localMomentumCosTheta = -1;
	}
	double localMomentumSinTheta = sqrt(1 - localMomentumCosTheta*localMomentumCosTheta);

	localMomentumX = localMomentum*localMomentumSinTheta*cos(localMomentumPhi);
	localMomentumY = localMomentum*localMomentumSinTheta*sin(localMomentumPhi);
	localMomentumZ = localMomentum*localMomentumCosTheta;

	isCosmicRay = false;
	setAbsoluteMomentum(bin);
	initialMomentum = absoluteMomentum;
	previousAbsoluteMomentum = initialMomentum;
	//path = std::list<double>();
	writePath = wpath;

}

Particle::~Particle(){
	//path.clear();
}

void Particle::setAbsoluteMomentum(double U, double Utheta, double Uphi){
	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = localMomentum*localMomentum;
	double V = sqrt(sqrp/(sqrm+sqrp/sqrc));
	double theta = acos(localMomentumZ/localMomentum);
	if(abs(localMomentum) < DBL_EPSILON){
		theta = pi/2;
	}
	double phi = atan2(localMomentumY,localMomentumX);

	double ux;
	double uy;
	double uz;

	if((thetagridNumber == 1)&&(phigridNumber == 1)){
		ux = 0;
		uy = 0;
		uz = U;		
	} else {

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

	if(abs(1 - absoluteV*absoluteV/c2) < DBL_EPSILON){
		printf("absluteV = c in setAbsoluteMmentum\n");
	}
	absoluteMomentum = mass*absoluteV/sqrt(1 - absoluteV*absoluteV/c2);
	if(absoluteMomentum != absoluteMomentum){
		printf("absoluteMomentum != absoluteMomentum\n");
	}

	absoluteMomentumTheta = acos(particleAbsoluteV.z/absoluteV);
	if(abs(absoluteV) < DBL_EPSILON){
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
		ux = 0;
		uy = 0;
		uz = U;		
		//if( U < 0){
		//	U = U;
		//}
	} else {
		//double c2 = speed_of_light*speed_of_light;

		ux = U*sin(Utheta)*cos(Uphi);
		uy = U*sin(Utheta)*sin(Uphi);
		uz = U*cos(Utheta);
	}

	//double localV1 = getLocalV();

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
			printf("v > c\n");
		}
	}

	if(abs(1 - localV*localV/c2) < DBL_EPSILON){
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
	if(abs(localMomentum) < DBL_EPSILON){
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

		ux = bin->U*sin(bin->UTheta)*cos(bin->UPhi);
		uy = bin->U*sin(bin->UTheta)*sin(bin->UPhi);
		uz = bin->U*cos(bin->UTheta);
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

	if(abs(1 - absoluteV*absoluteV/c2) < DBL_EPSILON){
		printf("absoluteV = c in setAbsoluteMomentum\n");
	}
	absoluteMomentum = mass*absoluteV/sqrt(1 - absoluteV*absoluteV/c2);
	if(absoluteMomentum != absoluteMomentum){
		printf("absoluteMomentum != absoluteMomentum\n");
	}

	absoluteMomentumTheta = acos(particleAbsoluteV.z/absoluteV);
	if(abs(absoluteV) < DBL_EPSILON){
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

		ux = bin->U*sin(bin->UTheta)*cos(bin->UPhi);
		uy = bin->U*sin(bin->UTheta)*sin(bin->UPhi);
		uz = bin->U*cos(bin->UTheta);
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


	if(abs(1 - localV*localV/c2) < DBL_EPSILON){
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
		printf("v != v\n");
	}
	return v;
}

double Particle::getLocalV(){
	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = localMomentum*localMomentum;
	double v = sqrt(sqrp/(sqrm+sqrp/sqrc));
	if( v != v){
		printf("v != v\n");
	}
	return v;
}

double Particle::getAbsoluteR(){
	return sqrt(absoluteX*absoluteX + absoluteY*absoluteY + absoluteZ*absoluteZ);
}

double Particle::getAbsoluteTheta(){
	double r = getAbsoluteR();
	if(abs(r) < DBL_EPSILON){
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
	/*double v = getAbsoluteV();
	//double c2 = speed_of_light*speed_of_light;
	if( (v*v)/c2 >= 1){
		printf(" v > c in getEnergy()\n");
	}*/
	return sqrt(mass*mass*c2*c2 + absoluteMomentum*absoluteMomentum*c2) - mass*c2;
}

double Particle::getInitialEnergy(){
	/*double v = getAbsoluteV();
	//double c2 = speed_of_light*speed_of_light;
	if( (v*v)/c2 >= 1){
		printf(" v > c in getEnergy()\n");
	}*/
	return sqrt(mass*mass*c2*c2 + initialMomentum*initialMomentum*c2) - mass*c2;
}

double Particle::getRadialSpeed(){
	return getAbsoluteV()*(cos(absoluteMomentumPhi)*sin(absoluteMomentumTheta)*cos(getAbsolutePhi())*sin(getAbsoluteTheta()) + sin(absoluteMomentumPhi)*sin(absoluteMomentumTheta)*sin(getAbsolutePhi())*sin(getAbsoluteTheta()) + cos(absoluteMomentumTheta)*cos(getAbsoluteTheta()));
}

void Particle::moveToBinRight(SpaceBin* bin){
	absoluteZ = bin->r1 + (bin->r2 - bin->r1)*(1-epsilon);
}