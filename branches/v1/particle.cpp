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
	path = std::list<double>();
	writePath = wPath;
}

Particle::Particle(int a, int znumber, SpaceBin* bin, bool wpath){
	double x = bin->r2*2*(uniRandom() - 0.5);
	double y = bin->r2*2*(uniRandom() - 0.5);
	double z = bin->r2*2*(uniRandom() - 0.5);
	double theta = acos(z/sqrt(x*x + y*y +z*z));
	double phi = atan2(y,x);
	if(phi < 0){
		phi = phi + 2*pi;
	}
	double r = sqrt(x*x + y*y + z*z);
	while(!(order(bin->r1,r,bin->r2) && order(bin->theta1, theta, bin->theta2) && order(bin->phi1, phi, bin->phi2))){
		x = bin->r2*2*(uniRandom() - 0.5);
		y = bin->r2*2*(uniRandom() - 0.5);
		z = bin->r2*2*(uniRandom() - 0.5);
		theta = acos(z/sqrt(x*x + y*y +z*z));
		phi = atan2(y,x);
		if(phi < 0){
			phi = phi + 2*pi;
		}
		r = sqrt(x*x + y*y + z*z);
	}

	absoluteX = x;
	absoluteY = y;
	absoluteZ = z;
	A = a;
	Z = znumber;
	mass = A*massProton;
	weight = bin->density*bin->volume/mass;

	double px = randomMaxwell(bin->temperature, mass);
	double py = randomMaxwell(bin->temperature, mass);
	double pz = randomMaxwell(bin->temperature, mass);

	localMomentum = sqrt(px*px+py*py+pz*pz);
	initialLocalMomentum = localMomentum;
	localMomentumX = px;
	localMomentumY = py;
	localMomentumZ = pz;
	isCosmicRay = false;
	setAbsoluteMomentum(bin->U,bin->UTheta,bin->UPhi);
	initialMomentum = absoluteMomentum;
	path = std::list<double>();
	writePath = wpath;

}

/*void Particle::LorentzTransition(double v1,double vtheta1, double vphi1, double v2,double vtheta2, double vphi2){

	double v1z = v1*cos(vtheta1);
	double v1x = v1*sin(vtheta1)*cos(vphi1);
	double v1y = v1*sin(vtheta1)*sin(vphi1);

	double v2z = v2*cos(vtheta1);
	double v2x = v2*sin(vtheta1)*cos(vphi1);
	double v2y = v2*sin(vtheta1)*sin(vphi1);

	Matrix3d matrix = Matrix3d::createBasisByOneVector(vector3d(v1x,v1y,v1z));

	vector3d v2rot = matrix.multiply(vector3d(v2x,v2y,v2z));


	vector3d u = summVelocity(v2rot,v1);

	Matrix3d invertMatrix = matrix.Inverse();

	u = invertMatrix.multiply(u);

	matrix = Matrix3d.createBasisByOneVector(u);

	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = localMomentum*localMomentum;

	double v = sqrt(sqrp/(sqrm+sqrp/sqrc));
	double sqrv = v*v;

	vector3d particleV = vector3d(v*sin(localMomentumTheta)*cos(localMomentumPhi),v*sin(localMomentumTheta)*sin(localMomentumPhi),v*cos(localMomentumTheta));

	particleV = matrix.multiply(particleV);



	double vx2 = (vx1-u)/(1-vx1*u/sqrc);
	double vy2 = vy1*sqrt(1-u*u/sqrc)/(1-vx1*u/sqrc);
	localPitchAngle = atan2(vy2,vx2);
	v = sqrt(vx2*vx2+vy2*vy2);
	localMomentum = mass*v/sqrt(1 - v*v/sqrc);
	//double localMomentumX = vx2*mass/sqrt(1-vx2*vx2/(sqrc*cos(localPitchAngle)*cos(localPitchAngle)));
	/////////////for debug;
	setAbsoluteMomentum(v2);
}*/

void Particle::setAbsoluteMomentum(double U, double Utheta, double Uphi){
	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = localMomentum*localMomentum;
	double V = sqrt(sqrp/(sqrm+sqrp/sqrc));
	//double c2 = speed_of_light*speed_of_light;
	double theta = acos(localMomentumZ/localMomentum);
	double phi = atan2(localMomentumY,localMomentumX);

	double ux = U*sin(Utheta)*cos(Uphi);
	double uy = U*sin(Utheta)*sin(Uphi);
	double uz = U*cos(Utheta);
	
	//u - скорость плазмы в абсолютной СО. Локальная ось z направлена по скорости плазмы
	Matrix3d* matrix = Matrix3d::createBasisByOneVector(vector3d(ux,uy,uz));

	vector3d particleV = vector3d(V*sin(theta)*cos(phi),V*sin(theta)*sin(phi),V*cos(theta));

	vector3d particleAbsoluteV = summVelocity(particleV,-U);

	Matrix3d* invertMatrix = matrix->Inverse();

	double absoluteV = particleAbsoluteV.getNorm();

	particleAbsoluteV = matrix->multiply(particleAbsoluteV);

	absoluteV = particleAbsoluteV.getNorm();

	if(absoluteV > speed_of_light){
		if(absoluteV < (1 + epsilon)*speed_of_light){
			absoluteV = (1 - epsilon)*speed_of_light;
			printf("v = c");
		} else {
			printf("aaa");
		}
	}

	absoluteMomentum = mass*absoluteV/sqrt(1 - absoluteV*absoluteV/c2);
	if(absoluteMomentum != absoluteMomentum){
		printf("aaa");
	}

	absoluteMomentumTheta = acos(particleAbsoluteV.z/absoluteV);
	absoluteMomentumPhi = atan2(particleAbsoluteV.y, particleAbsoluteV.x);
	
	delete matrix;
	delete invertMatrix;
}

void Particle::setLocalMomentum(double U, double Utheta,double Uphi){
	double sqrm = mass*mass;
	double sqrc = speed_of_light*speed_of_light;
	double sqrp = absoluteMomentum*absoluteMomentum;
	double V = sqrt(sqrp/(sqrm+sqrp/sqrc));
	//double c2 = speed_of_light*speed_of_light;
	double theta = absoluteMomentumTheta;
	double phi = absoluteMomentumPhi;

	double ux = U*sin(Utheta)*cos(Uphi);
	double uy = U*sin(Utheta)*sin(Uphi);
	double uz = U*cos(Utheta);

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
			printf("aaa");
		}
	}

	localMomentum = mass*localV/sqrt(1 - localV*localV/c2);

	//if(abs(localMomentum - initialLocalMomentum) > epsilon*localMomentum){
		//printf("aaa");
	//}

	if(localMomentum != localMomentum){
		printf("aaa");
	}
	localMomentumX = mass*particleLocalV.x/sqrt(1 - localV*localV/c2);
	localMomentumY = mass*particleLocalV.y/sqrt(1 - localV*localV/c2);
	localMomentumZ = mass*particleLocalV.z/sqrt(1 - localV*localV/c2);

	delete matrix;
	delete invertMatrix;
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
		printf("aaa");
	}
	return v;
}

double Particle::getAbsoluteR(){
	return sqrt(absoluteX*absoluteX + absoluteY*absoluteY + absoluteZ*absoluteZ);
}

double Particle::getAbsoluteTheta(){
	return acos(absoluteZ/getAbsoluteR());
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