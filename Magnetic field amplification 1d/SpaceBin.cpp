#include "stdafx.h"
#include "SpaceBin.h"
#include "particle.h"
#include "constants.h"
#include "simulation.h"
#include "util.h"
#include "BinFlux.h"


SpaceBin::SpaceBin(){
	magneticField = new double[kgridNumber];
	for (int i = 0; i < kgridNumber;++i){
		magneticField[i] = 0;
	}
	sortedParticles = new std::list<Particle*>[kgridNumber];
}
SpaceBin::SpaceBin(double R, double Theta, double Phi, double deltar, double deltatheta, double deltaphi, double u, double rho, double utheta,double uphi, double t, double b, int i, int j, int k){
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
	UTheta = utheta;
	UPhi = uphi;

	density = rho;
	temperature = t;

	B0 = b;

	numberR = i;
	numberTheta = j;
	numberPhi = k;

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
		printf("aaa");
	}
	particle->setLocalMomentum(U,UTheta,UPhi);
	double lambda = getFreePath(particle);
	if(lambda != lambda){
		printf("aaa");
	}
	if(lambda == 0){
		printf("bbbbbbb");
	}
	//double c2 = speed_of_light*speed_of_light;
	double gammaFactor = 1/sqrt(1 - U*U/c2);
	double colisionTime = lambda/particle->getLocalV();
	int l = 0;
	while(isInBin(particle) && (time < timeStep)){
		//падает где-то здесь
		l++;
		if(l > 1000){
			//printf("%d \n",l);
		}
	//while(((particle.absoluteR > r1)&&(particle.absoluteR < r2))&&((particle.absoluteTheta > theta1)&&(particle.absoluteTheta < theta2))&&((particle.absolutePhi > phi1)&&(particle.absolutePhi < phi2))){
		if(colisionTime > gammaFactor*(timeStep - time)){
			colisionTime = gammaFactor*(timeStep - time);
		}
		time += colisionTime/gammaFactor;
		if(time != time){
			printf("aaa");
		}
		if(colisionTime != colisionTime){
			printf("aaa");
		}
		makeOneStep(particle, colisionTime);
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
	double phi = atan2(particle->absoluteY, particle->absoluteX);
	if(phi < 0) {
		phi = phi + 2*pi;
	}


	return binByCoordinates(particle->absoluteZ, theta, phi,0, r2 - r1, theta2 - theta1, phi2 - phi1);
}

int* SpaceBin::binByCoordinates(double r, double theta, double phi, double r0, double deltar, double deltatheta, double deltaphi){
	int i = lowerInt((r - r0)/deltar);
	int j = lowerInt(theta/deltatheta);
	int k = lowerInt(phi/deltaphi);
	if(i > rgridNumber - 1){
		//printf("aaa");
	}
	if( i < 0){
		//printf("aaa");
	}
	if(j > thetagridNumber - 1){
		j = lowerInt((theta - epsilon)/deltatheta);
		if(j > thetagridNumber - 1){
			printf("aaa");
		}
	}
	if(k > phigridNumber - 1){
		printf("aaa");
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
		printf("aaa");
	}
	if( 0*lambda != 0*lambda){
		printf("aaa");
	}
	return lambda;
}

void SpaceBin::makeOneStep(Particle* particle, double colisionTime){
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

	double maxTheta = sqrt(6*deltat/colisionTime);

	particle->absoluteX = tempx;
	particle->absoluteY = tempy;
	particle->absoluteZ = tempz;

	scattering(particle, maxTheta);
}

void SpaceBin::scattering(Particle* particle, double maxTheta){
	double localTheta = acos(particle->localMomentumZ/particle->localMomentum);
	double localPhi = atan2(particle->localMomentumY, particle->localMomentumX);
	double deltaTheta = acos(1-uniRandom()*(1-cos(maxTheta)));
	double deltaPhi = uniRandom()*2*pi;
	localTheta = acos(cos(localTheta)*cos(deltaTheta)+sin(localTheta)*sin(deltaTheta)*cos(deltaPhi));
	double sinLocalTheta = sin(localTheta);
	localPhi = localPhi - asin(sin(deltaPhi)*sin(deltaTheta)/sinLocalTheta);
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
		//неправда какая-то
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
	for(int j = 0; j < kgridNumber; ++j){
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
	/*double u;
	double t;
	double rho;

	if (number < zeroPoint){
		//double discriminant = (crMomentaFlux - momentaFlux0)*(crMomentaFlux - momentaFlux0)-4*(energyFlux0 - crEnergyFlux)*(massFlux0 - crMassFlux)*(1-1/gamma)*(1-(gamma-1)/(2*gamma));
		//if (discriminant > 0){
			///u = (momentaFlux0 - crMomentaFlux + sqrt(discriminant))/(2*(massFlux0 - crMassFlux)*(1-(gamma-1)/(2*gamma)));
			//Вопрос + или -?
			//u = (momentaFlux0 - crMomentaFlux + sqrt(discriminant))/(2*(massFlux0 - crMassFlux)*(1-(gamma-1)/(2*gamma)));
			u = -(momentaFlux - momentaFlux0)/(density0*U0);
			//if (u < 0){
				//printf("u < 0");
			//}
			u = (1-weight)*U+weight*u;
			//t = (energyFlux0 - crEnergyFlux- (massFlux0 - crMassFlux)*u*u/2)/((massFlux0 - crMassFlux)*(1-1/gamma)*kBoltzman/massProton);
			rho = particleMassFlux/u;

		//U = (1-weight)*U+weight*u;
		//temperature = (1-weight)*temperature+weight*t;
		//density= (1-weight)*density+weight*rho;
			if ((u > 0)&&(u < speed_of_light)){
				if(abs(U-u) < 0.2*U){
					U=u;
				} else {
					if(u > U){
						U = 1.2*U;
					} else {
						U = 0.8*U;
					}
				}
			} else {
				U = 0.8*U;
			}
			//if (t > 0){
				//temperature = t;
			//}
			if (rho > 0){
				//TODO!
				//density = rho;
			}
		//} else {
			//printf("%s", "warning! discriminant < 0, do nothing\n");
		//}
	//} else {

	}*/
}
////интегрирует спектральную плотность и вычисляет эффективное поле 
void SpaceBin::updateMagneticField(){
	W =0;
	double dk = (maxK - minK)/(kgridNumber - 1);
	for(int i = 0; i < kgridNumber; ++i){
		W +=magneticField[i]*dk;
	}
	/////////// TODO 4 или 8?
	B =sqrt(B0*B0 + 4*pi*W);
}

//// Возвращает плотность давления, создаваемого частицами в близости резонанса
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
	detectedParticlesR1.push_back(new Particle(*particle));
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
	detectedParticlesR2.push_back(new Particle(*particle));
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

	


