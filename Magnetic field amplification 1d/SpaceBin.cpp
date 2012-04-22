#include "stdafx.h"
#include "SpaceBin.h"
#include "particle.h"
#include "constants.h"
#include "simulation.h"
#include "util.h"
#include "BinFlux.h"
#include "output.h"


SpaceBin::SpaceBin(){
	magneticField = new double[kgridNumber];
	for (int i = 0; i < kgridNumber;++i){
		magneticField[i] = 0;
	}
	sortedParticles = new std::list<Particle*>[kgridNumber];
	averageVelocity = 0;
	cosmicRayBoundMomentum = 0;
}
SpaceBin::SpaceBin(double R, double Theta, double Phi, double deltar, double deltatheta, double deltaphi, double u, double rho, double utheta,double uphi, double t, double b, int i, int j, int k, bool smallAngle){
	r = R;
	theta = Theta;
	phi = Phi;
	r1 = r - deltar/2;
	r2 = r + deltar/2;
	theta1 = theta - deltatheta/2;
	theta2 = theta + deltatheta/2;
	phi1 = phi - deltaphi/2;
	phi2 = phi + deltaphi/2;

	volume = r2 - r1;

	U = u;
	massVelocity = u;
	averageVelocity = u;
	UTheta = utheta;
	UPhi = uphi;

	density = rho;
	temperature = t;

	B0 = b;

	numberR = i;
	numberTheta = j;
	numberPhi = k;

	smallAngleScattering = smallAngle;

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

	particleMassFlux.reset();
	particleMomentaFlux.reset();
	particleEnergyFlux.reset();

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
	invertMatrix = matrix->Inverse();
	cosmicRayBoundMomentum = 0;
}

SpaceBin::~SpaceBin(){
	resetDetectors();
	delete sortedParticles;
	delete[] magneticField;
}

int* SpaceBin::propagateParticle(Particle* particle ,double& time, double timeStep){
	if(time != time){
		printf("time != time");
	}
	//particle->setLocalMomentum(U,UTheta,UPhi);
	particle->setLocalMomentum(this);
	double lambda = getFreePath(particle);
	if(lambda != lambda){
		printf("lambda != lambda");
	}
	if(lambda == 0){
		printf("lambda == 0");
	}
	//double c2 = speed_of_light*speed_of_light;
	if(abs(1 - U*U/c2) < DBL_EPSILON){
		printf("u = c in propageteParticle\n");
	}
	double gammaFactor = 1/sqrt(1 - U*U/c2);
	double colisionTime = getFreeTime(particle);
	if(smallAngleScattering){
	int l = 0;
	while(isInBin(particle) && (time < timeStep)){
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
			printf("time != time");
		}
		if(colisionTime != colisionTime){
			printf("colisionTime != colisionTime");
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

	double particleR = sqrt(particle->absoluteX*particle->absoluteX + particle->absoluteY*particle->absoluteY + particle->absoluteZ*particle->absoluteZ);
	double theta = acos(particle->absoluteZ/particleR);
	if( abs(r) < DBL_EPSILON){
		theta = pi/2;
	}
	double phi = atan2(particle->absoluteY, particle->absoluteX);
	if(phi < 0) {
		phi = phi + 2*pi;
	}

	//particle->setAbsoluteMomentum(U,UTheta,UPhi);
	particle->setAbsoluteMomentum(this);
	return binByCoordinates(particle->absoluteZ, theta, phi,0, r2 - r1, theta2 - theta1, phi2 - phi1);
	} else {
		//double v = particle->getLocalV()*particle->localMomentumZ/particle->localMomentum;
		//if(particle->localMomentum < epsilon*sqrt(particle->mass*temperature*kBoltzman)){
			double absoluteV = particle->getAbsoluteV()*cos(particle->absoluteMomentumTheta);
			double deltat;
			if(absoluteV > 0){
				deltat = (r2 - particle->absoluteZ)/absoluteV; 
				if(deltat < gammaFactor*(timeStep - time)){
					particle->absoluteZ = r2*(1 + epsilon);
					time += deltat/gammaFactor;
				} else {
					particle->absoluteZ += absoluteV*(timeStep - time)*gammaFactor;
					deltat = (timeStep - time)*gammaFactor;
					time = timeStep;
				}
			} else {
				deltat = (r1 - particle->absoluteZ)/absoluteV;
				if(deltat < gammaFactor*(timeStep - time)){
					particle->absoluteZ = r1*(1 - epsilon);
					if(r1 < epsilon){
						particle->absoluteZ = -epsilon;
					}
					time += deltat/gammaFactor;
				} else {
					particle->absoluteZ += absoluteV*(timeStep - time)*gammaFactor;
					deltat = (timeStep - time)*gammaFactor;
					time = timeStep;
				}
			}
			particle->absoluteX += particle->getAbsoluteV()*sin(particle->absoluteMomentumTheta)*cos(particle->absoluteMomentumPhi);
			particle->absoluteY += particle->getAbsoluteV()*sin(particle->absoluteMomentumTheta)*sin(particle->absoluteMomentumPhi);
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
			double particleR = sqrt(particle->absoluteX*particle->absoluteX + particle->absoluteY*particle->absoluteY + particle->absoluteZ*particle->absoluteZ);
			double theta = acos(particle->absoluteZ/particleR);
			if( abs(r) < DBL_EPSILON){
				theta = pi/2;
			}
			double phi = atan2(particle->absoluteY, particle->absoluteX);
			if(phi < 0) {
				phi = phi + 2*pi;
			}
			return binByCoordinates(particle->absoluteZ, theta, phi,0, r2 - r1, theta2 - theta1, phi2 - phi1);
		//}
	}
}

int* SpaceBin::binByCoordinates(double r, double theta, double phi, double r0, double deltar, double deltatheta, double deltaphi){
	int i = lowerInt((r - r0)/deltar);
	int j = lowerInt(theta/deltatheta);
	if( theta == pi ){
		theta = thetagridNumber - 1;
	}
	int k = lowerInt(phi/deltaphi);
	if (phi == 2*pi){
		k = 0;
	}
	if(i > rgridNumber - 1){
		//printf("aaa");
	}
	if( i < 0){
		//printf("aaa");
	}
	if(j > thetagridNumber - 1){
		j = lowerInt((theta - epsilon)/deltatheta);
		if(j > thetagridNumber - 1){
			printf("j > thetagridNumber - 1");
		}
	}
	if(k > phigridNumber - 1){
		printf("k > phigridNumber - 1");
		k = lowerInt(phi/deltaphi);
	}
	int* result = new int[3];
	result[0] = i;
	result[1] = j;
	result[2] = k;
	return result;
}

double SpaceBin::getFreePath(Particle* particle){
	//return speed_of_light*particle.localMomentum/(particle.Z*electron_charge*B);
	double lambda = speed_of_light*particle->localMomentum/(particle->Z*electron_charge*B0);
	if( lambda != lambda){
		printf("lambda != lambda");
	}
	if( 0*lambda != 0*lambda){
		printf("0*lambda != 0*lambda");
	}
	return lambda;
}

double SpaceBin::getFreeTime(Particle* particle){
	double lambda = getFreePath(particle);
	if(abs(particle->getLocalV()) > DBL_EPSILON){
		return lambda/particle->getLocalV();
	} else {
		printf("particle localv = 0 in getFreeTime\n");
		return speed_of_light*particle->mass/(particle->Z*electron_charge*B0);
	}
}

void SpaceBin::makeOneStep(Particle* particle, double colisionTime, double& time){
	double deltat = colisionTime/defaultTimeRelation;
	double particleU = particle->getAbsoluteV();
	double localMomentum = particle->localMomentum;
	double localMomentumZ = particle->localMomentumZ;

	double r = particle->absoluteZ;
	double theta = particle->getAbsoluteTheta();
	double phi = particle->getAbsolutePhi();

	double ux = particleU*sin(particle->absoluteMomentumTheta)*cos(particle->absoluteMomentumPhi);
	double uy = particleU*sin(particle->absoluteMomentumTheta)*sin(particle->absoluteMomentumPhi);
	double uz = particleU*cos(particle->absoluteMomentumTheta);

	double tempx = particle->absoluteX + ux*deltat;
	double tempy = particle->absoluteY + uy*deltat;
	double tempz = particle->absoluteZ + uz*deltat;

	double tempr = tempz;
	//double temptheta = acos(tempz/tempr);
	//double tempphi = atan2(tempy, tempx);
	//if(tempphi < 0){
		//tempphi = tempphi + 2*pi;
	//}

	/*double Ur = (sin(particle->getAbsoluteTheta())*sin(particle->absoluteMomentumTheta)*cos(particle->getAbsolutePhi())*cos(particle->absoluteMomentumPhi) + 
		sin(particle->getAbsoluteTheta())*sin(particle->absoluteMomentumTheta)*sin(particle->getAbsolutePhi())*sin(particle->absoluteMomentumPhi) +
		cos(particle->absoluteMomentumTheta)*cos(particle->getAbsoluteTheta()))*particleU;*/
	/*double Utheta = ( cos(particle.absoluteTheta)*sin(particle.absoluteMomentumTheta)*cos(particle.absolutePhi)*cos(particle.absoluteMomentumPhi) +
		cos(particle.absoluteTheta)*sin(particle.absoluteMomentumTheta)*sin(particle.absolutePhi)*sin(particle.absoluteMomentumPhi) -
		sin(particle.absoluteTheta)*cos(particle.absoluteMomentumTheta))*particleU;
	double Uphi = (- sin(particle.absolutePhi)*sin(particle.absoluteMomentumTheta)*cos(particle.absoluteMomentumPhi) +
		cos(particle.absolutePhi)*sin(particle.absoluteMomentumTheta)*sin(particle.absoluteMomentumPhi))*particleU;*/
	//TODO tempphi > 2*pi!!!!!
	//double deltaPhi =  angleDelta(phi, tempphi);
	//double deltaTheta =  temptheta - theta;
	double deltaR =  tempr - r;

	//if(deltaR < 0){
		//printf("bbb");
	//}

	//if(deltaPhi != deltaPhi){
		//printf("aaaa");
	//}
	//if(0*deltaPhi != 0*deltaPhi){
		//printf("aaaa");
	//}

	while (abs(deltaR) > (r2 - r1)/10) {
		deltat = deltat/10;
		tempx = particle->absoluteX + ux*deltat;
		tempy = particle->absoluteY + uy*deltat;
		tempz = particle->absoluteZ + uz*deltat;

		tempr = tempz;

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
	if(particle->localMomentum == 0.0){
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
	if(abs(sinDeltaPhi) > 1){
		printf("%s %lf\n","sinDeltaPhi > 1",sinDeltaPhi);
		localPhi = 0;
	}
	if(abs(sinLocalTheta) < DBL_EPSILON){
		localPhi = 0;
	}
	particle->localMomentumZ = particle->localMomentum*cos(localTheta);
	particle->localMomentumX = particle->localMomentum*sinLocalTheta*cos(localPhi);
	particle->localMomentumY = particle->localMomentum*sinLocalTheta*sin(localPhi);
	particle->setAbsoluteMomentum(this);
}

void SpaceBin::updateFluxes(){
	//updateCosmicRayFluxes();
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
			if(abs(1 - (v*v)/c2) < DBL_EPSILON){
				printf("v = c in updateCosmicRayFluxes\n");
			}
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
			if(abs(1 - (v*v)/c2) < DBL_EPSILON){
				printf("v = c in updateCosmicRayFluxes\n");
			}
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
	std::list<Particle*>::iterator it = detectedParticlesR2.begin();
	while(it != detectedParticlesR2.end()){
		Particle* particle = *it;
		delete particle;
		++it;
	}
	it = detectedParticlesR1.begin();
	while(it != detectedParticlesR1.end()){
		Particle* particle = *it;
		delete particle;
		++it;
	}
	/*it = particles.begin();
	while(it != particles.end()){
		Particle* particle = *it;
		delete particle;
		++it;
	}*/
	particles.clear();
	detectedParticlesR1.clear();
	detectedParticlesR2.clear();
	detectedParticlesTheta1.clear();
	detectedParticlesTheta2.clear();
	detectedParticlesPhi1.clear();
	detectedParticlesPhi2.clear();
	for(int j = 0; j < kgridNumber; ++j){
		it = sortedParticles[j].begin();
		while(it != sortedParticles[j].end()){
			Particle* particle = *it;
			delete particle;
			++it;
		}
		sortedParticles[j].clear();
	}
	massFlux.fluxR1 = 0;
	massFlux.fluxR2 = 0;
	momentaFlux.fluxR1 = 0;
	momentaFlux.fluxR2 = 0;
	energyFlux.fluxR1 = 0;
	energyFlux.fluxR2 = 0;
	particleMassFlux.fluxR1 = 0;
	particleMassFlux.fluxR2 = 0;
	particleMomentaFlux.fluxR1 = 0;
	particleMomentaFlux.fluxR2 = 0;
	particleEnergyFlux.fluxR1 = 0;
	particleEnergyFlux.fluxR2 = 0;
	crMassFlux.fluxR1 = 0;
	crMassFlux.fluxR2 = 0;

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
	double r = particle->absoluteZ;
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
	if(numberR == febNumber){
		//detectedParticlesR1.push_back(new Particle(*particle));
	}
	//double c2 = speed_of_light*speed_of_light;
	if(particle->isCosmicRay){
		crMassFlux.fluxR1 += particle->mass*particle->weight;
	}
	particleMassFlux.fluxR1 += particle->mass*particle->weight;
	double v = particle->getAbsoluteV();
	double vr = v*cos(particle->absoluteMomentumTheta); 
	particleMomentaFlux.fluxR1 += vr*particle->mass*particle->weight/sqrt(1 - (v*v)/c2);
	if(abs(1 - (v*v)/c2) < DBL_EPSILON){
		printf("v = c in detectParticleR1\n");
	}
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
	double vr = v*cos(particle->absoluteMomentumTheta);
	if(abs(1 - (v*v)/c2) < DBL_EPSILON){
		printf("v = c in detectParticleR2\n");
	}
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

void SpaceBin::resetVelocity(std::list<Particle*>& particles){
	//if(averageVelocity < 0){
	//	averageVelocity = averageVelocity;
	//}
	resetVelocity(particles, -2*abs(averageVelocity), 2*abs(averageVelocity));
}

void SpaceBin::resetVelocity(std::list<Particle*>& particles, double u1, double u2){
	if(abs(u1-u2) < epsilon*abs(u2)){
		U = (u1 + u2)/2;
		return;
	}
	double u = (u1 + u2)/2;

	if(u >= 0){
		if( evaluateMomentumFlux(particles, u) > 0){
			resetVelocity(particles,u,u2);
		} else {
			resetVelocity(particles, u1, u);
		}
	} else {
		if( evaluateMomentumFlux(particles, u) > 0){
			resetVelocity(particles,u1,u);
		} else {
			resetVelocity(particles, u, u2);
		}
	}
}

double SpaceBin::evaluateMomentumFlux(std::list<Particle*>& particles, double u){
	double result = 0;
	std::list<Particle*>::iterator it = particles.begin();
	while(it != particles.end()){
		Particle* particle = *it;
		particle->setLocalMomentum(u,0,0);
		result += particle->localMomentumZ*particle->weight/getFreeTime(particle);
		++it;
	}
	return result;
}

void SpaceBin::updateCosmicRayBoundMomentum(){

	/*std::list<Particle*>::iterator it = particles.begin();

	if(it == particles.end()){
		return;
	}

	double tempMomentum;
	double* averageVz = new double[pgridNumber];
	double* count = new double[pgridNumber];
	
	for(int i = 0; i < pgridNumber; ++i){
		averageVz[i] = 0.0;
		count[i] = 0.0;
	}

	Particle* particle = *it;

	double minp = particle->localMomentum;
	double maxp = minp;

	while( it != particles.end()){
		particle = *it;
		if(maxp < particle->localMomentum){
			maxp = particle->localMomentum;
		}
		if(minp > particle->localMomentum){
			minp = particle->localMomentum;
		}
		++it;
	}

	double deltap = (maxp - minp)/(pgridNumber );

	it = particles.begin();
	while( it != particles.end()){
		Particle* particle = *it;
		double p = particle->localMomentum;
		int i = lowerInt((p - minp)/deltap);
		if(( i < 0) || (i >= pgridNumber)){
			//printf("index out of bound in updateCosmicRayBoundMomentum\n");
			if(i >= pgridNumber){
				i = pgridNumber - 1;
			}
		}
		if(particle->localMomentum > DBL_EPSILON*sqrt(massProton*kBoltzman*temperature)){
			averageVz[i] += particle->getLocalV()*particle->weight*particle->localMomentumZ/particle->localMomentum;
			count[i] += particle->weight;
		} else {
			printf("localMomentum = 0\n");
		}
		++it;
	}

	tempMomentum = 0;
	double Vz = 0;
	for(int i = pgridNumber - 1; i >=0; --i){
		Vz += averageVz[i];
		if(Vz > 0){
		//if(averageVz[i] > 0){
			tempMomentum = minp + (i + 1)*deltap;
			break;
		}
	}

	for(int i = 0; i < pgridNumber; ++i){
		if(count[i] > 0){
			averageVz[i] /= count[i];
		}
	}
	
	cosmicRayBoundMomentum = tempMomentum;

	double j  = 0;

	it = particles.begin();
	while(it != particles.end()){
		Particle* particle = *it;
		if(particle->localMomentum > cosmicRayBoundMomentum){
			j += electron_charge*particle->Z*(particle->weight/volume)*particle->getLocalV()*particle->localMomentumZ/(particle->localMomentum*particle->mass);
		}
		++it;
	}

	crFlux = j;

	//const char* name = "./output/averageVz" + numberR;

	if(numberR == 10){
		outputAverageVz(minp, maxp, averageVz, "./output/averageVz1.dat");
	}
	if(numberR == 20){
		outputAverageVz(minp, maxp, averageVz, "./output/averageVz2.dat");
	}
	if(numberR == 30){
		outputAverageVz(minp, maxp, averageVz, "./output/averageVz3.dat");
	}
	if(numberR == 40){
		outputAverageVz(minp, maxp, averageVz, "./output/averageVz4.dat");
	}
	if(numberR == 50){
		outputAverageVz(minp, maxp, averageVz, "./output/averageVz5.dat");
	}
	if(numberR == 60){
		outputAverageVz(minp, maxp, averageVz, "./output/averageVz6.dat");
	}
	if(numberR == 70){
		outputAverageVz(minp, maxp, averageVz, "./output/averageVz7.dat");
	}
	if(numberR == 80){
		outputAverageVz(minp, maxp, averageVz, "./output/averageVz8.dat");
	}
	if(numberR == 90){
		outputAverageVz(minp, maxp, averageVz, "./output/averageVz9.dat");
	}
	if(numberR == 0){
		outputAverageVz(minp, maxp, averageVz, "./output/averageVz0.dat");
	}
	delete[] averageVz;
	delete[] count;*/
	if( particles.size() > 0){
		std::list<Particle*>::iterator it = particles.begin();
		double maxp;
		maxp = (*it)->localMomentumZ;
		++it;
		double weight = 0;
		weight += (*it)->weight;
		while(it != particles.end()){
			Particle* particle = *it;
			weight += particle->weight;
			if(maxp < abs(particle->localMomentumZ)){
				maxp = abs(particle->localMomentumZ);
			}
			++it;
		}
		double* distribution = new double[pgridNumber];
		for (int i = 0; i < pgridNumber; ++i){
			distribution[i] = 0;
		}
		int particleNumber = 0;
		double mass;
		double deltap = (2*maxp)/(pgridNumber);
		it = particles.begin();
		while(it != particles.end()){
			particleNumber++;
			Particle* particle = *it;
			mass = particle->mass;
			++it;
			double p = particle->localMomentumZ;
			if (abs(p) >= maxp){
				continue;
			}
			for(int i =0; i< pgridNumber; ++i){
				if (p <-maxp + deltap*(i + 1)){
					distribution[i] += particle->weight/(weight*deltap);
					break;
				}
			}
		}


		/*double sigma;
		int zeroP = pgridNumber/2;
		int sigmaIndex = zeroP;
		double maxDistribution = (distribution[zeroP] + distribution[zeroP + 1])/2;
		for(int i = zeroP;i < pgridNumber-2; ++i){
			if((distribution[i - 1] + distribution[i] + distribution[i + 1])/3< maxDistribution/exp(0.5)){
				sigmaIndex = i - 1;
				break;
			}
		}
		double d1 = distribution[sigmaIndex];
		double d2 = distribution[sigmaIndex + 1];

		sigma =  (maxDistribution/exp(0.5) - d1)*deltap/(d2 - d1) + (sigmaIndex + 0.5 - zeroP)*deltap;
		temperature = sigma*sigma/(massProton*kBoltzman);*/


		updateTemperature(distribution, deltap);

		if(numberR == 10){
			FILE* outPDF = fopen("./output/zpdf1.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(-maxp + i*deltap), distribution[i], maxwell(-maxp + (i + 1/2)*deltap, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 20){
			FILE* outPDF = fopen("./output/zpdf2.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(-maxp + i*deltap), distribution[i], maxwell(-maxp + (i + 1/2)*deltap, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 30){
			FILE* outPDF = fopen("./output/zpdf3.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(-maxp + i*deltap), distribution[i], maxwell(-maxp + (i + 1/2)*deltap, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 40){
			FILE* outPDF = fopen("./output/zpdf4.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(-maxp + i*deltap), distribution[i], maxwell(-maxp + (i + 1/2)*deltap, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 50){
			FILE* outPDF = fopen("./output/zpdf5.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(-maxp + i*deltap), distribution[i], maxwell(-maxp + (i + 1/2)*deltap, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 60){
			FILE* outPDF = fopen("./output/zpdf6.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(-maxp + i*deltap), distribution[i], maxwell(-maxp + (i + 1/2)*deltap, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 70){
			FILE* outPDF = fopen("./output/zpdf7.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(-maxp + i*deltap), distribution[i], maxwell(-maxp + (i + 1/2)*deltap, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 80){
			FILE* outPDF = fopen("./output/zpdf8.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(-maxp + i*deltap), distribution[i], maxwell(-maxp + (i + 1/2)*deltap, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 90){
			FILE* outPDF = fopen("./output/zpdf9.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(-maxp + i*deltap), distribution[i], maxwell(-maxp + (i + 1/2)*deltap, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 0){
			FILE* outPDF = fopen("./output/zpdf0.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(-maxp + i*deltap), distribution[i], maxwell(-maxp + (i + 1/2)*deltap, massProton, temperature));
			}
			fclose(outPDF);
		}

		/*FILE* outPDF = fopen(fileName,"w");
		for(int i = 0; i < pgridNumber; ++i){
			fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(-maxp + i*deltap), distribution[i], maxwell(-maxp + (i + 1/2)*deltap, massProton, temperature));
		}
		fclose(outPDF);*/


		double j = 0;

		for(int i = 0; i < pgridNumber; ++i){
			distribution[i] -= maxwell(-maxp + (i + 1/2)*deltap, massProton, temperature);
			j += weight*((-maxp + (i + 1/2)*deltap)/massProton*volume)*electron_charge*distribution[i];
		}

		crFlux = j;

		delete[] distribution;
	}
}

void SpaceBin::updateTemperature(double* distribution, double deltap){
	double sigma;
	int zeroP = pgridNumber/2;
	int sigmaIndex = zeroP;
	double maxDistribution = (distribution[zeroP] + distribution[zeroP + 1])/2;
	for(int i = zeroP;i < pgridNumber-2; ++i){
		if((distribution[i - 1] + distribution[i] + distribution[i + 1])/3< maxDistribution/exp(0.5)){
			sigmaIndex = i - 1;
			break;
		}
	}
	double d1 = distribution[sigmaIndex];
	double d2 = distribution[sigmaIndex + 1];


	sigma =  (maxDistribution/exp(0.5) - d1)*deltap/(d2 - d1) + (sigmaIndex + 0.5 - zeroP)*deltap;
	temperature = sigma*sigma/(massProton*kBoltzman);

	double minT = 0;
	double maxT = (deltap*pgridNumber)*(deltap*pgridNumber)/(8*kBoltzman*massProton);
	double T = findTemperature(minT,maxT,distribution, deltap);
	temperature = T;
}

double SpaceBin::findTemperature(double minT, double maxT, double* distribution, double deltap){
	double tempT = (minT + maxT)/2;
	double T1 = (minT + tempT)/2;
	double T2 = (tempT + maxT)/2;

	if((T1 - T2) < 100){
		return tempT;
	}

	double delta1 = getError(distribution, deltap, T1);
	double delta2 = getError(distribution, deltap, T2);

	if(delta1 < delta2){
		return findTemperature(minT, tempT, distribution, deltap);
	} else {
		return findTemperature(tempT,maxT, distribution, deltap);
	}
}

double SpaceBin::getError(double* distribution, double deltap, double T){
	double error = 0;
	for(int i = 0; i < pgridNumber/2 - 2; ++i){
		double p = i*deltap;
		error += abs((distribution[pgridNumber/2 + i - 1] + distribution[pgridNumber/2 + i] + distribution[pgridNumber/2 + i + 1])/3 - maxwell((i + 1/2)*deltap, massProton, T));
	}
	return error;
}

	


