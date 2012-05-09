#include "stdafx.h"
#include "SpaceBin.h"
#include "particle.h"
#include "constants.h"
#include "simulation.h"
#include "util.h"
#include "BinFlux.h"
#include "matrix3d.h"


SpaceBin::SpaceBin(){
	magneticField = new double[kgridNumber];
	for (int i = 0; i < kgridNumber;++i){
		magneticField[i] = 0;
	}
	sortedParticles = new std::list<Particle*>[kgridNumber];
}
SpaceBin::SpaceBin(double R, double Theta, double Phi, double deltar, double deltatheta, double deltaphi, double u, double rho, double utheta,double uphi, double t, double b, int i, int j, int k, bool scattering){
	r = R;
	theta = Theta;
	phi = Phi;
	r1 = r - deltar/2;
	r2 = r + deltar/2;
	theta1 = theta - deltatheta/2;
	theta2 = theta + deltatheta/2;
	phi1 = phi - deltaphi/2;
	phi2 = phi + deltaphi/2;

	volume = (r2*r2*r2 - r1*r1*r1)*(phi2 - phi1)*(cos(theta1) - cos(theta2))/3;

	U = u;
	UTheta = utheta;
	UPhi = uphi;
	averageVelocity = U;

	density = rho;
	temperature = t;

	B0 = b;

	numberR = i;
	numberTheta = j;
	numberPhi = k;

	smallAngleScattering = scattering;

	double ux;
	double uy;
	double uz;
	if((thetagridNumber == 1)&&(phigridNumber == 1)){
		ux = 0;
		uy = 0;
		uz = U;		
	} else {
		//double c2 = speed_of_light*speed_of_light;

		ux = U*sin(UTheta)*cos(UPhi);
		uy = U*sin(UTheta)*sin(UPhi);
		uz = U*cos(UTheta);
	}	
	
	//u - скорость плазмы в абсолютной —ќ. Ћокальна€ ось z направлена по скорости плазмы
	matrix = Matrix3d::createBasisByOneVector(vector3d(ux,uy,uz));
	if((0*matrix->matrix[0][0] != 0*matrix->matrix[0][0])
		||(0*matrix->matrix[0][1] != 0*matrix->matrix[0][1])
		||(0*matrix->matrix[0][2] != 0*matrix->matrix[0][2])
		||(0*matrix->matrix[1][0] != 0*matrix->matrix[1][0])
		||(0*matrix->matrix[1][1] != 0*matrix->matrix[1][1])
		||(0*matrix->matrix[1][2] != 0*matrix->matrix[1][2])
		||(0*matrix->matrix[2][0] != 0*matrix->matrix[2][0])
		||(0*matrix->matrix[2][1] != 0*matrix->matrix[2][1])
		||(0*matrix->matrix[2][2] != 0*matrix->matrix[2][2])){
		printf("NaN matrix\n");
	}
	invertMatrix = matrix->Inverse();

	//TODO!!!!!!!!!!!!!!!!
	/*magneticField = new double[kgridNumber];
	for (int i = 0; i < kgridNumber;++i){
		magneticField[i] = 0;
	}*/
	detectedParticlesR1 = std::list<Particle*>();
	detectedParticlesR2 = std::list<Particle*>();
	detectedParticlesTheta1 = std::list<Particle*>();
	detectedParticlesTheta2 = std::list<Particle*>();
	detectedParticlesPhi1 = std::list<Particle*>();
	detectedParticlesPhi2 = std::list<Particle*>();
	sortedParticles = new std::list<Particle*>[kgridNumber];
}

SpaceBin::~SpaceBin(){
	delete[] magneticField;
}

int* SpaceBin::propagateParticle(Particle* particle ,double& time, double timeStep){
	if(time != time){
		printf("aaa\n");
	}
	particle->setLocalMomentum(this);
	double lambda = getFreePath(particle);
	if(lambda != lambda){
		printf("aaa\n");
	}
	if(lambda == 0){
		printf("lambda == 0\n");
	}
	//double c2 = speed_of_light*speed_of_light;
	double gammaFactor = 1/sqrt(1 - U*U/c2);
	double colisionTime = lambda/particle->getLocalV();
	int l = 0;
	if(smallAngleScattering){
		while(isInBin(particle) && (time < timeStep)){
			//падает где-то здесь
			l++;
			if(l > 1000){
				//printf("%d \n",l);
			}
		//while(((particle.absoluteR > r1)&&(particle.absoluteR < r2))&&((particle.absoluteTheta > theta1)&&(particle.absoluteTheta < theta2))&&((particle.absolutePhi > phi1)&&(particle.absolutePhi < phi2))){
			if(colisionTime > defaultTimeRelation*gammaFactor*(timeStep - time)){
				colisionTime = defaultTimeRelation*gammaFactor*(timeStep - time);
			}
			//time += colisionTime/(defaultTimeRelation*gammaFactor);
			if(time != time){
				printf("aaa\n");
			}
			if(colisionTime != colisionTime){
				printf("aaa\n");
			}
			makeOneStep(particle, colisionTime, time);
		}
		/*if(particle.absoluteTheta < 0){
			particle.absoluteTheta = -particle.absoluteTheta;
			particle.absolutePhi += pi;
		} 
		if(particle.absoluteTheta > pi){
			particle.absoluteTheta = 2*pi - particle.absoluteTheta;
			particle.absolutePhi += pi;
		}
		while(particle.absolutePhi < 0){
			particle.absolutePhi += 2*pi;
		}
		while(particle.absolutePhi > 2*pi){
			particle.absolutePhi -= 2*pi;
		}*/

		double r = sqrt(particle->absoluteX*particle->absoluteX + particle->absoluteY*particle->absoluteY + particle->absoluteZ*particle->absoluteZ);
		double theta = acos(particle->absoluteZ/r);
		if( abs(r) < epsilon){
			theta = pi/2;
		}
		double phi = atan2(particle->absoluteY, particle->absoluteX);
		if(phi < 0) {
			phi = phi + 2*pi;
		}


		return binByCoordinates(r, theta, phi,r1 - numberR*(r2 - r1), r2 - r1, theta2 - theta1, phi2 - phi1);
	} else {
		largeAngleScattering(particle, time, timeStep);
		double particleR = sqrt(particle->absoluteX*particle->absoluteX + particle->absoluteY*particle->absoluteY + particle->absoluteZ*particle->absoluteZ);
		double theta = acos(particle->absoluteZ/particleR);
		if( abs(r) < epsilon){
			theta = pi/2;
		}
		double phi = atan2(particle->absoluteY, particle->absoluteX);
		if(phi < 0) {
			phi = phi + 2*pi;
		}
		return binByCoordinates(particleR, theta, phi,0, r2 - r1, theta2 - theta1, phi2 - phi1);
	}
}

void SpaceBin::largeAngleScattering(Particle* particle, double& time, double timeStep){
	double gammaFactor = 1/sqrt(1 - U*U/c2);
	double lambda = getFreePath(particle);
	double colisionTime = lambda/particle->getLocalV();
	double localV = particle->getLocalV();
	if(localV < epsilon){
		colisionTime = timeStep*2;
	}
	double r0;
	double sign;
	double absoluteV = particle->getAbsoluteV();

	double b = particle->absoluteX*absoluteV*sin(particle->absoluteMomentumTheta)*cos(particle->absoluteMomentumPhi)+ particle->absoluteY*absoluteV*sin(particle->absoluteMomentumTheta)*sin(particle->absoluteMomentumPhi)  + particle->absoluteZ*absoluteV*cos(particle->absoluteMomentumTheta);

	double discr1 = sqr(b) - sqr(absoluteV)*(sqr(particle->getAbsoluteR()) - r1*r1);
	double discr2 = sqr(b) - sqr(absoluteV)*(sqr(particle->getAbsoluteR()) - r2*r2);
	double discr;
	if((b < 0) && (discr1 > 0) && (numberR > 0)){
		discr = discr1;
		sign = -1;
	} else {
		discr = discr2;
		sign = 1;
	}
	if(discr < 0){
		printf("discr < 0\n");
	}

	if(absoluteV < epsilon){
		printf("absoluteV = 0 in largeAngleScattering\n");
	}

	double deltat1 = (-b - sqrt(discr))/(absoluteV*absoluteV);
	double deltat2 = (-b + sqrt(discr))/(absoluteV*absoluteV);

	double deltat;
	if(deltat1 <= 0){
		deltat = deltat2;
	} else {
		deltat = deltat1;
	}

	if(deltat > gammaFactor*(timeStep - time)){
		deltat = gammaFactor*(timeStep - time);
	}
	if(deltat != deltat){
		printf("deltat != deltat in largeAngleScattering\n");
	}
	if(0*deltat != 0*deltat){
		printf("deltat = infinity in largeAngleScattering\n");
	}
	time += deltat/gammaFactor;

	if(deltat < 0){
		printf("dt < 0\n");
		FILE* errFile = fopen("./output/errFile.dat","a");
		fprintf(errFile,"dt < 0\n");
		fprintf(errFile,"x= %lf, y= %lf, z= %lf, vx= %lf, vy= %lf, vz= %lf, r1= %lf, r2=%lf\n",particle->absoluteX,particle->absoluteY,particle->absoluteZ,absoluteV*sin(particle->absoluteMomentumTheta)*cos(particle->absoluteMomentumPhi),absoluteV*sin(particle->absoluteMomentumTheta)*sin(particle->absoluteMomentumPhi),absoluteV*cos(particle->absoluteMomentumTheta),r1,r2);
		fclose(errFile);
	}

	if(numberR == 0){
		sign = 1;
	}

	particle->absoluteX += (absoluteV*sin(particle->absoluteMomentumTheta)*cos(particle->absoluteMomentumPhi)*deltat)*(1+sign*epsilon);
	particle->absoluteY += absoluteV*sin(particle->absoluteMomentumTheta)*sin(particle->absoluteMomentumPhi)*deltat*(1+sign*epsilon);
	particle->absoluteZ += absoluteV*cos(particle->absoluteMomentumTheta)*deltat*(1+sign*epsilon);

	double probability;
	if(deltat > colisionTime){
		probability = 1;
	} else {
		probability = deltat/colisionTime;
	}
	if(uniRandom() <= probability){
		double cosLocalTheta = 2*(uniRandom() - 0.5);
		double phi = 2*pi*uniRandom();
		double theta = acos(cosLocalTheta);
		particle->localMomentumZ = particle->localMomentum*cos(theta);
		particle->localMomentumX = particle->localMomentum*sin(theta)*cos(phi);
		particle->localMomentumY = particle->localMomentum*sin(theta)*sin(phi);
		particle->setAbsoluteMomentum(this);
	}
}


int* SpaceBin::binByCoordinates(double r, double theta, double phi, double r0, double deltar, double deltatheta, double deltaphi){
	int i = lowerInt((r - r0)/deltar);
	int j = lowerInt(theta/deltatheta);
	if( theta == pi ){
		theta = thetagridNumber - 1;
	}
	int k = lowerInt(phi/deltaphi);
	if(k >= phigridNumber){
		k = 0;
	}
	if(k < 0){
		k = phigridNumber - 1;
	}
	if (phi == 2*pi){
		k = 0;
	}
	if(i > rgridNumber - 1){
		//printf("aaa");
	}
	if( i < 0){
		//printf("aaa");
	}
	if(j >= thetagridNumber){
		j = 0;
	}
	if(j < 0){
		j = thetagridNumber - 1;
	}
	if(k > phigridNumber - 1){
		printf("aaa");
	}
	int* result = new int[3];
	result[0] = i;
	//result[1] = j;
	//result[2] = k;
	//TODO !!!!!!!!!
	result[1]=j;
	result[2]=k;
	return result;
}

double SpaceBin::getFreePath(Particle* particle){
	//return speed_of_light*particle.localMomentum/(particle.Z*electron_charge*B);
	double lambda = speed_of_light*particle->localMomentum/(particle->Z*electron_charge*B0);
	if( lambda != lambda){
		printf("aaa");
	}
	if( 0*lambda != 0*lambda){
		printf("aaa");
	}
	return lambda;
}

void SpaceBin::makeOneStep(Particle* particle, double colisionTime, double& time){
	double deltat = colisionTime/defaultTimeRelation;
	double particleU = particle->getAbsoluteV();
	double localMomentum = particle->localMomentum;
	double localMomentumZ = particle->localMomentumZ;

	double r = particle->getAbsoluteR();
	double theta = particle->getAbsoluteTheta();
	double phi = particle->getAbsolutePhi();

	double ux = particleU*sin(particle->absoluteMomentumTheta)*cos(particle->absoluteMomentumPhi);
	double uy = particleU*sin(particle->absoluteMomentumTheta)*sin(particle->absoluteMomentumPhi);
	double uz = particleU*cos(particle->absoluteMomentumTheta);

	//deltat = min4(deltat, abs(particle->absoluteX/ux), abs(particle->absoluteY/uy), abs(particle->absoluteZ/uz));
	if(deltat < epsilon){
		deltat = epsilon;
	}

	double tempx = particle->absoluteX + ux*deltat;
	double tempy = particle->absoluteY + uy*deltat;
	double tempz = particle->absoluteZ + uz*deltat;

	double tempr = sqrt(tempx*tempx + tempy*tempy +tempz*tempz);
	double temptheta = acos(tempz/tempr);
	if(abs(tempr) < epsilon){
		temptheta = theta;
	}
	double tempphi = atan2(tempy, tempx);
	if(tempphi < 0.0){
		tempphi = tempphi + 2*pi;
	}

	/*double Ur = (sin(particle->getAbsoluteTheta())*sin(particle->absoluteMomentumTheta)*cos(particle->getAbsolutePhi())*cos(particle->absoluteMomentumPhi) + 
		sin(particle->getAbsoluteTheta())*sin(particle->absoluteMomentumTheta)*sin(particle->getAbsolutePhi())*sin(particle->absoluteMomentumPhi) +
		cos(particle->absoluteMomentumTheta)*cos(particle->getAbsoluteTheta()))*particleU;*/
	/*double Utheta = ( cos(particle.absoluteTheta)*sin(particle.absoluteMomentumTheta)*cos(particle.absolutePhi)*cos(particle.absoluteMomentumPhi) +
		cos(particle.absoluteTheta)*sin(particle.absoluteMomentumTheta)*sin(particle.absolutePhi)*sin(particle.absoluteMomentumPhi) -
		sin(particle.absoluteTheta)*cos(particle.absoluteMomentumTheta))*particleU;
	double Uphi = (- sin(particle.absolutePhi)*sin(particle.absoluteMomentumTheta)*cos(particle.absoluteMomentumPhi) +
		cos(particle.absolutePhi)*sin(particle.absoluteMomentumTheta)*sin(particle.absoluteMomentumPhi))*particleU;*/
	//TODO tempphi > 2*pi!!!!!
	double deltaPhi =  angleDelta(phi, tempphi);
	double deltaTheta =  temptheta - theta;
	double deltaR =  tempr - r;

	if(deltaR < 0){
		//printf("bbb");
	}

	if(deltaPhi != deltaPhi){
		printf("aaaa");
	}
	if(0*deltaPhi != 0*deltaPhi){
		printf("aaaa");
	}

	while( ((phigridNumber > 1)&&(abs(deltaPhi) > pi/(10*phigridNumber))) || ((thetagridNumber > 1)&&(abs(deltaTheta) > pi/(2*10*thetagridNumber))) || (abs(deltaR) > (r2 - r1)/10) || (!(isInThisOrNear(tempr, temptheta, tempphi))) || ( temptheta > pi) || (temptheta < 0)){
		deltat = deltat/10;
		tempx = particle->absoluteX + ux*deltat;
		tempy = particle->absoluteY + uy*deltat;
		tempz = particle->absoluteZ + uz*deltat;

		tempr = sqrt(tempx*tempx + tempy*tempy +tempz*tempz);
		temptheta = acos(tempz/tempr);
		if(abs(tempr) < epsilon){
			temptheta = theta;
		}
		tempphi = atan2(tempy, tempx);
		if(tempphi < 0){
			tempphi = tempphi + 2*pi;
		}
		deltaPhi =  angleDelta(phi, tempphi);
		deltaTheta =  temptheta - theta;
		deltaR =  tempr - r;

		if(deltaR < 0){
			//printf("bbb");
		}
	}

	time += deltat;

	double maxTheta = sqrt(6*deltat/colisionTime);

	particle->absoluteX = tempx;
	particle->absoluteY = tempy;
	particle->absoluteZ = tempz;

	scattering(particle, maxTheta);
}

void SpaceBin::scattering(Particle* particle, double maxTheta){
	double localTheta = acos(particle->localMomentumZ/particle->localMomentum);
	if(abs(particle->localMomentum) < epsilon){
		localTheta = pi/2;
	}
	double localPhi = atan2(particle->localMomentumY, particle->localMomentumX);
	double cosDeltaTheta = 1-uniRandom()*(1-cos(maxTheta));
	if(abs(cosDeltaTheta > 1)){
		printf("cosDeltaTheta > 1\n");
	}
	double deltaTheta = acos(cosDeltaTheta);
	double deltaPhi = uniRandom()*2*pi;
	double cosLocalTheta = cos(localTheta)*cosDeltaTheta+sin(localTheta)*sin(deltaTheta)*cos(deltaPhi);
	localTheta = acos(cosLocalTheta);
	if( cosLocalTheta > 1){
		printf("cosLocalTheta > 1\n");
		localTheta = 0;
	}
	if( cosLocalTheta < -1){
		printf("cosLocalTheta < -1\n");
		localTheta = pi;
	}
	double sinLocalTheta = sin(localTheta);
	double sinDeltaPhi = sin(deltaPhi)*sin(deltaTheta)/sinLocalTheta;
	localPhi = localPhi - asin(sinDeltaPhi);
	if(abs(sinDeltaPhi > 1)){
		printf("%s %lf\n","sinDeltaPhi > 1",sinDeltaPhi);
		localPhi = 0;
	}
	if(abs(sinLocalTheta) < epsilon){
		localPhi = 0;
	}
	particle->localMomentumZ = particle->localMomentum*cos(localTheta);
	particle->localMomentumX = particle->localMomentum*sinLocalTheta*cos(localPhi);
	particle->localMomentumY = particle->localMomentum*sinLocalTheta*sin(localPhi);
	particle->setAbsoluteMomentum(U,UTheta,UPhi);
}

void SpaceBin::updateFluxes(){
	updateCosmicRayFluxes();
	/*updateBulkFluxes();
	updateThermalFluxes();
	updateMagneticFluxes();


	massFlux = particleMassFlux+bulkMassFlux;
	momentaFlux = particleMomentaFlux+bulkMomentaFlux+thMomentaFlux+magMomentaflux;
	energyFlux = particleEnergyFlux+bulkEnergyFlux+thEnergyFlux+magEnergyFlux;
	if (momentaFlux > 1000000000000000000000.0){
		printf("aaaa");
	}*/
}
////////// I don't sure how to do it
void SpaceBin::updateCosmicRayFluxes(){
	//particleMassFlux.reset();
	//particleMomentaFlux.reset();
	//particleEnergyFlux.reset();
	//crMassFlux.reset();
	//double c2 = speed_of_light*speed_of_light;
	if((thetagridNumber == 1) && (phigridNumber == 1)){
		std::list<Particle*>::iterator it = detectedParticlesR1.begin();
		/*if(number < zeroPoint){
			density = 0;
		}*/
		while(it != detectedParticlesR1.end()){
			Particle particle = **it;
			//crMassFlux += density0*partOfCosmicRay*particle.getAbsoluteVx(U);
			//double n;
			//double rho;
			//double v = particle.getAbsoluteVx();
			if(particle.isCosmicRay){
				crMassFlux.fluxR1 += particle.mass*particle.weight;
			}
			particleMassFlux.fluxR1 += particle.mass*particle.weight;
			double v = particle.getAbsoluteV();
			double vr = particle.getRadialSpeed(); 
			particleMomentaFlux.fluxR1 += vr*particle.mass*particle.weight/sqrt(1 - (v*v)/c2);
			particleEnergyFlux.fluxR1 += particle.weight*particle.getEnergy();
		}
		it = detectedParticlesR2.begin();
		while(it != detectedParticlesR2.end()){
			Particle particle = **it;
			if(particle.isCosmicRay){
				crMassFlux.fluxR2 += particle.mass*particle.weight;
			}
			particleMassFlux.fluxR2 += particle.mass*particle.weight;
			double v = particle.getAbsoluteV();
			double vr = particle.getRadialSpeed(); 
			particleMomentaFlux.fluxR2 += vr*vr*particle.mass*particle.weight/sqrt(1 - (v*v)/c2);
			particleEnergyFlux.fluxR2 += particle.weight*particle.getEnergy();
		}
		//if(v < 0){
		//	printf("v < 0");
		//}
		//неправда кака€-то
		//n = (particle.weight*density0/(massProton*particle.A))*abs(particle.initialSpeedX/v);
		//n = abs(density0*upstreamSpeed*particle.weight/v);
		//if(number < zeroPoint){
			//density += n*massProton*particle.A;
		//}
		//particleMassFlux += n*massProton*particle.A*v;
		//particleMomentaFlux += n*particle.getMomentumX()*v;
		//particleEnergyFlux += n*c2*(1/sqrt(1-v*v/c2)-1)*particle.A*massProton*v;
		//if( particle.isCosmicRay){
			//crMassFlux += n*massProton*particle.A*v;
		//}
		//crMomentaFlux += (density0*partOfCosmicRay/(massProton*particle.A)) * particle.momentumX*particle.getAbsoluteVx(U);
		//crEnergyFlux += (density0*partOfCosmicRay)* c2*(1/sqrt(1-v*v/c2)-1) *particle.getAbsoluteVx(U);
		//++it;
	}
	//double momentaflux0 = density*U*U;

}

void SpaceBin::updateBulkFluxes(){
	/*bulkMassFlux = 0;
	bulkMomentaFlux = 0;
	bulkEnergyFlux = 0;
	if(number > zeroPoint){
		bulkMassFlux = density*U;
		bulkMomentaFlux = bulkMassFlux*U;
		bulkEnergyFlux = bulkMomentaFlux*U/2;
	}
	//bulkMassFlux = 0;
	//bulkMomentaFlux = 0;
	//bulkEnergyFlux = 0;*/
}
void SpaceBin::updateThermalFluxes(){
	//thMomentaFlux = density/massProton*kBoltzman*temperature;
	//thEnergyFlux = thMomentaFlux*gamma/(gamma-1)*U;
	//thMomentaFlux = 0;
	//thEnergyFlux = 0;
}
void SpaceBin::resetDetectors(){
	detectedParticlesR1.clear();
	detectedParticlesR2.clear();
	detectedParticlesTheta1.clear();
	detectedParticlesTheta2.clear();
	detectedParticlesPhi1.clear();
	detectedParticlesPhi2.clear();

	std::list<Particle*>::iterator it = particles.begin();
	while(it != particles.end()){
		Particle* particle = *it;
		delete particle;
		++it;
	}
	particles.clear();

	for(int j = 0; j < kgridNumber; ++j){
		sortedParticles[j].clear();
	}
	massFlux.fluxR1 = 0;
	massFlux.fluxR2 = 0;
	massFlux.fluxTheta1 = 0;
	massFlux.fluxTheta2 = 0;
	massFlux.fluxPhi1 = 0;
	massFlux.fluxPhi2 = 0;
	momentaFlux.fluxR1 = 0;
	momentaFlux.fluxR2 = 0;
	momentaFlux.fluxTheta1 = 0;
	momentaFlux.fluxTheta2 = 0;
	momentaFlux.fluxPhi1 = 0;
	momentaFlux.fluxPhi2 = 0;
	energyFlux.fluxR1 = 0;
	energyFlux.fluxR2 = 0;
	energyFlux.fluxTheta1 = 0;
	energyFlux.fluxTheta2 = 0;
	energyFlux.fluxPhi1 = 0;
	energyFlux.fluxPhi2 = 0;
	particleMassFlux.fluxR1 = 0;
	particleMassFlux.fluxR2 = 0;
	particleMassFlux.fluxTheta1 = 0;
	particleMassFlux.fluxTheta2 = 0;
	particleMassFlux.fluxPhi1 = 0;
	particleMassFlux.fluxPhi2 = 0;
	particleMomentaFlux.fluxR1 = 0;
	particleMomentaFlux.fluxR2 = 0;
	particleMomentaFlux.fluxTheta1 = 0;
	particleMomentaFlux.fluxTheta2 = 0;
	particleMomentaFlux.fluxPhi1 = 0;
	particleMomentaFlux.fluxPhi2 = 0;
	particleEnergyFlux.fluxR1 = 0;
	particleEnergyFlux.fluxR2 = 0;
	particleEnergyFlux.fluxTheta1 = 0;
	particleEnergyFlux.fluxTheta2 = 0;
	particleEnergyFlux.fluxPhi1 = 0;
	particleEnergyFlux.fluxPhi2 = 0;
	crMassFlux.fluxR1 = 0;
	crMassFlux.fluxR2 = 0;
	crMassFlux.fluxTheta1 = 0;
	crMassFlux.fluxTheta2 = 0;
	crMassFlux.fluxPhi1 = 0;
	crMassFlux.fluxPhi2 = 0;

	double ux;
	double uy;
	double uz;
	if((thetagridNumber == 1)&&(phigridNumber == 1)){
		ux = 0;
		uy = 0;
		uz = U;		
	} else {
		//double c2 = speed_of_light*speed_of_light;

		ux = U*sin(UTheta)*cos(UPhi);
		uy = U*sin(UTheta)*sin(UPhi);
		uz = U*cos(UTheta);
	}	

	delete matrix;
	delete invertMatrix;
	matrix = Matrix3d::createBasisByOneVector(vector3d(ux,uy,uz));
	if((0*matrix->matrix[0][0] != 0*matrix->matrix[0][0])
		||(0*matrix->matrix[0][1] != 0*matrix->matrix[0][1])
		||(0*matrix->matrix[0][2] != 0*matrix->matrix[0][2])
		||(0*matrix->matrix[1][0] != 0*matrix->matrix[1][0])
		||(0*matrix->matrix[1][1] != 0*matrix->matrix[1][1])
		||(0*matrix->matrix[1][2] != 0*matrix->matrix[1][2])
		||(0*matrix->matrix[2][0] != 0*matrix->matrix[2][0])
		||(0*matrix->matrix[2][1] != 0*matrix->matrix[2][1])
		||(0*matrix->matrix[2][2] != 0*matrix->matrix[2][2])){
		printf("NaN matrix\n");
	}
	invertMatrix = matrix->Inverse();
	//momentaFlux.reset();
	//energyFlux.reset();
}

/////// I don't know how to do it
void SpaceBin::updateMagneticFluxes(){
	/*magEnergyFlux = (3/2)*U*W;
	////TODO как?
	magMomentaflux = W/2;
	//magEnergyFlux = 0;
	//magMomentaflux = 0;*/
}

void SpaceBin::multiplyParticleWeights(double value){
	/*std::list<Particle*>::iterator it = detectedParticlesR1.begin();
	std::list<Particle> new_list = std::list<Particle>();
	while(it != detectedParticlesR1.end()){
		Particle particle = *it;
		particle.weight *= value;
		new_list.push_back(particle);
		++it;
	}
	detectedParticlesR1 = new_list;

	it = detectedParticlesR2.begin();
	new_list = std::list<Particle>();
	while(it != detectedParticlesR2.end()){
		if(numberR == rgridNumber - 1){
			printf("aaa\n");
		}
		Particle particle = *it;
		particle.weight *= value;
		new_list.push_back(particle);
		++it;
	}
	detectedParticlesR2 = new_list;

	it = detectedParticlesTheta1.begin();
	new_list = std::list<Particle>();
	while(it != detectedParticlesTheta1.end()){
		Particle particle = *it;
		particle.weight *= value;
		new_list.push_back(particle);
		++it;
	}
	detectedParticlesTheta1 = new_list;

	it = detectedParticlesTheta2.begin();
	new_list = std::list<Particle>();
	while(it != detectedParticlesTheta2.end()){
		Particle particle = *it;
		particle.weight *= value;
		new_list.push_back(particle);
		++it;
	}
	detectedParticlesTheta2 = new_list;

	it = detectedParticlesPhi1.begin();
	new_list = std::list<Particle>();
	while(it != detectedParticlesPhi1.end()){
		Particle particle = *it;
		particle.weight *= value;
		new_list.push_back(particle);
		++it;
	}
	detectedParticlesPhi1 = new_list;

	it = detectedParticlesPhi2.begin();
	new_list = std::list<Particle>();
	while(it != detectedParticlesPhi2.end()){
		Particle particle = *it;
		particle.weight *= value;
		new_list.push_back(particle);
		++it;
	}
	detectedParticlesPhi2 = new_list;*/
}

void SpaceBin::resetProfile(double massFlux0, double momentaFlux0, double energyFlux0, double density0, double U0){

}
////интегрирует спектральную плотность и вычисл€ет эффективное поле 
void SpaceBin::updateMagneticField(){
	W =0;
	double dk = (maxK - minK)/(kgridNumber - 1);
	for(int i = 0; i < kgridNumber; ++i){
		W +=magneticField[i]*dk;
	}
	/////////// TODO 4 или 8?
	B =sqrt(B0*B0 + 4*pi*W);
}

//// ¬озвращает плотность давлени€, создаваемого частицами в близости резонанса
double SpaceBin::pressureSpectralDensity(int j, double density0, double U0, int Z, int A){
	/*double p = 0;
	double n = 0; 
	double rho = 0;
	double flux = 0;
	double k = minK + ((maxK - minK)*j)/(kgridNumber-1);
	//double pres = Z*electron_charge*B0/(speed_of_light*k);
	double kstep = (maxK - minK)/(kgridNumber-1);
	std::list<Particle>::iterator it = sortedParticles[j].begin();
	while ( it != sortedParticles[j].end()){
		Particle particle = *it;
		if (abs(particle.Z*electron_charge*B0/(speed_of_light*particle.localMomentum) - k) < kstep/2){
		//if (abs(particle.localMomentum-pres) < Z*electron_charge*B0*kstep/(2*speed_of_light*k*k)){
			double v = particle.getAbsoluteVx();
			double n1 = (particle.weight*density0/(massProton*particle.A))*abs(particle.initialSpeedX/v);
			n += n1;
			flux += n1*particle.getMomentumX()*v;
		}
		++it;
	}
	if(number < zeroPoint){
		rho = n*A*massProton;
		//TODO как считать p?
		p = flux;
		//p = flux - rho*U*U;
	} else {
		/////TODO ?
		//rho = xbins[i]->density;
		p = flux;
	}
	return p/kstep;*/
	return 0;
}

double SpaceBin::getMaxP(){
	std::list<Particle*>::iterator it = detectedParticlesR1.begin();
	double pmax = 0;
	while(it != detectedParticlesR1.end()){
		Particle particle = **it;
		if(particle.absoluteMomentum > pmax){
			pmax = particle.absoluteMomentum;
		}
		++it;
	}

	it = detectedParticlesR2.begin();
	while(it != detectedParticlesR2.end()){
		Particle particle = **it;
		if(particle.absoluteMomentum > pmax){
			pmax = particle.absoluteMomentum;
		}
		++it;
	}

	it = detectedParticlesTheta1.begin();
	while(it != detectedParticlesTheta1.end()){
		Particle particle = **it;
		if(particle.absoluteMomentum > pmax){
			pmax = particle.absoluteMomentum;
		}
		++it;
	}

	it = detectedParticlesTheta2.begin();
	while(it != detectedParticlesTheta2.end()){
		Particle particle = **it;
		if(particle.absoluteMomentum > pmax){
			pmax = particle.absoluteMomentum;
		}
		++it;
	}

	it = detectedParticlesPhi1.begin();
	while(it != detectedParticlesPhi1.end()){
		Particle particle = **it;
		if(particle.absoluteMomentum > pmax){
			pmax = particle.absoluteMomentum;
		}
		++it;
	}

	it = detectedParticlesPhi2.begin();
	while(it != detectedParticlesPhi2.end()){
		Particle particle = **it;
		if(particle.absoluteMomentum > pmax){
			pmax = particle.absoluteMomentum;
		}
		++it;
	}

	return pmax;
} 

double SpaceBin::getMinP(){
	double pmin = 1;
	std::list<Particle*>::iterator it = detectedParticlesR1.begin();
	if (it != detectedParticlesR1.end()){
		Particle particle = **it;
		pmin = particle.absoluteMomentum;
		++it;
		while(it != detectedParticlesR1.end()){
			Particle particle = **it;
			if(particle.absoluteMomentum < pmin){
				pmin = particle.absoluteMomentum;
			}
			++it;
		}
	}

	it = detectedParticlesR2.begin();
	while(it != detectedParticlesR2.end()){
		Particle particle = **it;
		if(particle.absoluteMomentum < pmin){
			pmin = particle.absoluteMomentum;
		}
		++it;
	}

	it = detectedParticlesTheta1.begin();
	while(it != detectedParticlesTheta1.end()){
		Particle particle = **it;
		if(particle.absoluteMomentum < pmin){
			pmin = particle.absoluteMomentum;
		}
		++it;
	}

	it = detectedParticlesTheta2.begin();
	while(it != detectedParticlesTheta2.end()){
		Particle particle = **it;
		if(particle.absoluteMomentum < pmin){
			pmin = particle.absoluteMomentum;
		}
		++it;
	}

	it = detectedParticlesPhi1.begin();
	while(it != detectedParticlesPhi1.end()){
		Particle particle = **it;
		if(particle.absoluteMomentum < pmin){
			pmin = particle.absoluteMomentum;
		}
		++it;
	}

	it = detectedParticlesPhi2.begin();
	while(it != detectedParticlesPhi2.end()){
		Particle particle = **it;
		if(particle.absoluteMomentum < pmin){
			pmin = particle.absoluteMomentum;
		}
		++it;
	}

	return pmin;
} 
double SpaceBin::getMaxLocalP(){
	std::list<Particle*>::iterator it = detectedParticlesR1.begin();
	double pmax = 0;
	while(it != detectedParticlesR1.end()){
		Particle particle = **it;
		if(particle.localMomentum > pmax){
			pmax = particle.localMomentum;
		}
		++it;
	}

	it = detectedParticlesR2.begin();
	while(it != detectedParticlesR2.end()){
		Particle particle = **it;
		if(particle.localMomentum > pmax){
			pmax = particle.localMomentum;
		}
		++it;
	}

	it = detectedParticlesTheta1.begin();
	while(it != detectedParticlesTheta1.end()){
		Particle particle = **it;
		if(particle.localMomentum > pmax){
			pmax = particle.localMomentum;
		}
		++it;
	}

	it = detectedParticlesTheta2.begin();
	while(it != detectedParticlesTheta2.end()){
		Particle particle = **it;
		if(particle.localMomentum > pmax){
			pmax = particle.localMomentum;
		}
		++it;
	}

	it = detectedParticlesPhi1.begin();
	while(it != detectedParticlesPhi1.end()){
		Particle particle = **it;
		if(particle.localMomentum > pmax){
			pmax = particle.localMomentum;
		}
		++it;
	}

	it = detectedParticlesPhi2.begin();
	while(it != detectedParticlesPhi2.end()){
		Particle particle = **it;
		if(particle.localMomentum > pmax){
			pmax = particle.localMomentum;
		}
		++it;
	}

	return pmax;
} 

double SpaceBin::getMinLocalP(){
	std::list<Particle*>::iterator it = detectedParticlesR1.begin();
	double pmin = 0;
	if (it != detectedParticlesR1.end()){
		Particle particle = **it;
		pmin = particle.localMomentum;
		++it;
		while(it != detectedParticlesR1.end()){
			Particle particle = **it;
			if(particle.localMomentum < pmin){
				pmin = particle.localMomentum;
			}
			++it;
		}
	}

	it = detectedParticlesR2.begin();
	while(it != detectedParticlesR2.end()){
		Particle particle = **it;
		if(particle.localMomentum < pmin){
			pmin = particle.localMomentum;
		}
		++it;
	}

	it = detectedParticlesTheta1.begin();
	while(it != detectedParticlesTheta1.end()){
		Particle particle = **it;
		if(particle.localMomentum < pmin){
			pmin = particle.localMomentum;
		}
		++it;
	}

	it = detectedParticlesTheta2.begin();
	while(it != detectedParticlesTheta2.end()){
		Particle particle = **it;
		if(particle.localMomentum < pmin){
			pmin = particle.localMomentum;
		}
		++it;
	}

	it = detectedParticlesPhi1.begin();
	while(it != detectedParticlesPhi1.end()){
		Particle particle = **it;
		if(particle.localMomentum < pmin){
			pmin = particle.localMomentum;
		}
		++it;
	}

	it = detectedParticlesPhi2.begin();
	while(it != detectedParticlesPhi2.end()){
		Particle particle = **it;
		if(particle.localMomentum < pmin){
			pmin = particle.localMomentum;
		}
		++it;
	}

	return pmin;
}

void SpaceBin::sortParticles(double minK, double maxK){
	double kstep = (maxK - minK)/(kgridNumber - 1);
	std::list<Particle*>::iterator it = detectedParticlesR2.begin();
	while (it != detectedParticlesR2.end()){
		int i = (int) ( (**it).localMomentum - minK)/kstep;
		sortedParticles[i].push_back(*it);
		++it;
	}
}

bool SpaceBin::isInBin(Particle* particle){
	double r = particle->getAbsoluteR();
	double theta = particle->getAbsoluteTheta();
	double phi = particle->getAbsolutePhi();
	if(r > r2) {
		//detectedParticlesR2.push_back(new Particle(*particle));
		detectParticleR2(particle);
		return false;
	}
	if(r < r1) {
		//detectedParticlesR1.push_back(new Particle(*particle));
		detectParticleR1(particle);
		return false;
	}
	if(thetagridNumber > 1){
		if(theta > theta2) {
			//detectedParticlesTheta2.push_back(new Particle(*particle));
			detectParticleTheta2(particle);
			return false;
		}
		if(theta < theta1){
			//detectedParticlesTheta1.push_back(new Particle(*particle));
			detectParticleTheta1(particle);
			return false;
		}
	}
	if(phigridNumber > 1){
		if(phi > phi2) {
			//detectedParticlesPhi2.push_back(new Particle(*particle));
			detectParticlePhi2(particle);
			return false;
		}
		if(phi < phi1){
			//detectedParticlesPhi1.push_back(new Particle(*particle));
			detectParticlePhi1(particle);
			return false;
		}
	}
	return true;
}

bool SpaceBin::isInThisOrNear(double r, double theta, double phi){
	if((thetagridNumber == 1) && (phigridNumber == 1)){
		return true;
	}
	if(thetagridNumber == 1) {
		return (order(r1, r, r2) || order(phi1, phi, phi2));
	}
	if(phigridNumber == 1){
		return (order(r1,r,r2) || order(theta1, theta, theta2));
	}
	if(order(r1, r, r2)){
		return (order(theta1, theta, theta2) || order(phi1, phi, phi2));
	} else {
		return (order(theta1, theta, theta2) && order(phi1, phi, phi2));
	}
}

void SpaceBin::detectParticleR1(Particle* particle){
	//detectedParticlesR1.push_back(new Particle(*particle));
	//double c2 = speed_of_light*speed_of_light;
	if(particle->isCosmicRay){
		crMassFlux.fluxR1 += particle->mass*particle->weight;
	}
	particleMassFlux.fluxR1 += particle->mass*particle->weight;
	double v = particle->getAbsoluteV();
	double vr = particle->getRadialSpeed(); 
	particleMomentaFlux.fluxR1 += vr*particle->mass*particle->weight/sqrt(1 - (v*v)/c2);
	particleEnergyFlux.fluxR1 += particle->weight*particle->getEnergy();
}

void SpaceBin::detectParticleR2(Particle* particle){
	//detectedParticlesR2.push_back(new Particle(*particle));
	//double c2 = speed_of_light*speed_of_light;
	if(particle->isCosmicRay){
		crMassFlux.fluxR2 += particle->mass*particle->weight;
	}
	particleMassFlux.fluxR2 += particle->mass*particle->weight;
	double v = particle->getAbsoluteV();
	double vr = particle->getRadialSpeed(); 
	particleMomentaFlux.fluxR2 += vr*particle->mass*particle->weight/sqrt(1 - (v*v)/c2);
	particleEnergyFlux.fluxR2 += particle->weight*particle->getEnergy();
}

void SpaceBin::detectParticleTheta1(Particle* particle){
}

void SpaceBin::detectParticleTheta2(Particle* particle){
}

void SpaceBin::detectParticlePhi1(Particle* particle){
}

void SpaceBin::detectParticlePhi2(Particle* particle){
}

	


