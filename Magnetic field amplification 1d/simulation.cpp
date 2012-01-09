#include "stdafx.h"
#include "simulation.h"
#include "SpaceBin.h"
#include "particle.h"
#include "output.h"
#include "util.h"
#include <omp.h>

Simulation::Simulation(){
	partOfCosmicRay = 0;
	allParticlesNumber = 0;
	kolmogorovCascading = true;
	resonantInstability = true;
	bellInstability = true;
	alpha = 1;
	beta = 1;
	gamma = 1;
	delta = 1;
	epsilon = 0.1;
	A = 1;
	Z = 1;
	startPDF = std::list <Particle*>();
	timeStep = defaultTimeStep;
	zeroBin = NULL;
	averageVelocity = new double[rgridNumber];
}

Simulation::~Simulation(){
	kolmogorovCascading = true;
	resonantInstability = true;
	bellInstability = true;
	alpha = 1;
	beta = 1;
	gamma = 1;
	delta = 1;
	epsilon = 0.1;
	delete[] averageVelocity;
}

void Simulation::initializeProfile(){
 	minK = defaultMinK;
	maxK = defaultMaxK;
	deltaR = (downstreamR - upstreamR)/(rgridNumber );
	deltaTheta = pi/(thetagridNumber );
	deltaPhi = 2*pi/(phigridNumber );
	double R = upstreamR + deltaR/2;
	double Theta = deltaTheta/2;
	double Phi = deltaPhi/2;
	bins = new SpaceBin***[rgridNumber];
	zeroBin = new SpaceBin(R - deltaR,Theta,Phi,deltaR,deltaTheta,deltaPhi,U0,density0,Theta,Phi,temperature,B0,-1,0,0);
    for(int i = 0; i < rgridNumber; ++i){
		bins[i] = new SpaceBin**[thetagridNumber];
		Theta = deltaTheta/2;
		for(int j = 0; j < thetagridNumber; ++j){
			bins[i][j] = new SpaceBin*[phigridNumber];
			Phi = deltaPhi/2;
			for(int k = 0; k < phigridNumber; ++k){
				double density = density0;
				double u;
				if(i < shockWavePoint){
					u = U0;
				} else {
					u = U0;
				}
				bins[i][j][k] = new SpaceBin(R,Theta,Phi,deltaR,deltaTheta,deltaPhi,u,density,Theta,Phi,temperature,B0,i,j,k);
				Phi = Phi + deltaPhi;
			}
			Theta = Theta + deltaTheta;
		}
		averageVelocity[i] = U0;
		R = R + deltaR;
	}
}
///// главная функция программы.
void Simulation::simulate(){
	FILE* outIteration;
	srand ( time(NULL) );
	initializeProfile();
	introducedParticles = getParticles();
	FILE* radialFile;
	outIteration = fopen("./output/tamc_iteration.dat","w");
	radialFile = fopen("./output/tamc_radial_profile.dat","w");
	fclose(outIteration);
	fclose(radialFile);
	for (int itNumber  = 0; itNumber < iterationNumber; ++itNumber){ 
		printf("%s", "\n");
		int j = 0;
		int l = 0;
		printf("%s", "Iteration started\n");
		if(itNumber == 0){
			printf("%s","First iteration");
			prevPoint = downstreamR;
			shockWavePoints.push_back(prevPoint);
			shockWaveVelocity.push_back(0.0);
		} else {
			printf("%s", "Particle propagation\n");
			std::list<Particle*>::iterator it = introducedParticles.begin();
			while( it != introducedParticles.end()){
				Particle* particle = *it;
				++it;
				//printf("%d %s",l,"\n");
				++l;
				bool side = false;
				double r = particle->getAbsoluteR();
				double theta = acos(particle->absoluteZ/r);
				double phi =atan2(particle->absoluteY, particle->absoluteX);
				if(phi < 0){
					phi = phi + 2*pi;
				}
				int* index = SpaceBin::binByCoordinates(particle->absoluteZ, theta, phi,upstreamR,deltaR,deltaTheta,deltaPhi);
				double time = 0;
				while((index[0] >= 0)){
					j++;
					SpaceBin* bin = bins[index[0]][index[1]][index[2]];
					if(time != time){
						printf("aaa");
					}
					int* tempIndex = bin->propagateParticle(particle,time, timeStep);
					delete[] index;
					index = tempIndex;
					if(index[0] > rgridNumber){
						printf("aaa");
					}
					if(index[0] == rgridNumber){
						//bin = bins[rgridNumber - 1][index[1]][index[2]];
						index[0] = rgridNumber - 1;
						bin->detectParticleR2(particle);
						particle->absoluteMomentumTheta = pi - particle->absoluteMomentumTheta;
						particle->absoluteMomentumPhi = 2*pi - particle->absoluteMomentumPhi;
						if(particle->absoluteMomentumPhi > 2*pi){
							particle->absoluteMomentumPhi -= 2*pi;
						}
						particle->moveToBinRight(bin);
						//bin->detectParticleR2(particle);
					}
					if (time >= timeStep){
						break;
					}
				}
				delete[] index;
			}
			/*printf("%s", "Fluxes updating\n");
			for (int i = 0; i < rgridNumber; ++i){
				for (int j = 0; j < phigridNumber; ++j){
					for (int k = 0; k < thetagridNumber; ++k){
						bins[i][j][k]->updateFluxes();
					}
				}
			}*/
			printf("%s", "Introducing new particles from left boundary\n");
			introduceNewParticles();

			printf("%s", "Reseting profile\n");
			resetProfile();
			collectAverageVelocity();
			resetVelocity();
			findShockWavePoint();
			printf("%s", "magnetic Field updating\n");
			//updateMagneticField();
			//output(*this);
			removeEscapedParticles();
			resetDetectors();
			printf("%s","iteration № ");
			printf("%d\n",itNumber);
		}
		outIteration = fopen("./output/tamc_iteration.dat","a");
		radialFile = fopen("./output/tamc_radial_profile.dat","a");
		fprintf(outIteration,"%s %d \n","iteration number ",itNumber);
		fclose(outIteration);
		double maxp;
		double minp;
		updateMaxMinP(minp, maxp);
		std::list<Particle*> list = std::list<Particle*>();
		for(int j = 0; j < thetagridNumber; ++j){
			for(int k = 0; k < phigridNumber; ++k){
				list.insert(list.end(),bins[rgridNumber - 1][j][k]->detectedParticlesR2.begin(),bins[rgridNumber - 1][j][k]->detectedParticlesR2.end());
			}
		}
		//outputPDF(list,"./output/tamc_pdf_down.dat",*this,minp,maxp);
		//outputStartPDF(list,"./output/tamc_pdf_start.dat",*this,minp,maxp);
		outputRadialProfile(bins,0,0,radialFile,averageVelocity);
		outputShockWave(shockWavePoints, shockWaveVelocity);
		fclose(radialFile);
	}
}

void Simulation::resetProfile(){
	if((thetagridNumber == 1) && (phigridNumber == 1)){
		//TODO what do with i = 0
		for(int i = 0; i < rgridNumber; i++){
			double massFlux1 = 0;
			double momentaFlux1 = 0;
			if(i > 0){
				massFlux1 = bins[i - 1][0][0]->particleMassFlux.fluxR2;
				momentaFlux1 = bins[i - 1][0][0]->particleMomentaFlux.fluxR2;
			}
			if(i == 0){
				massFlux1 = zeroBin->particleMassFlux.fluxR2;
				momentaFlux1 = zeroBin->particleMomentaFlux.fluxR2;
			}
			double massFlux2 = 0;
			double momentaFlux2 = 0;
			if ( i < rgridNumber - 1){
				massFlux2 = bins[i + 1][0][0]->particleMassFlux.fluxR1;
				momentaFlux2 = bins[i + 1][0][0]->particleMomentaFlux.fluxR1;
			}
			if(i == rgridNumber - 1){
				bins[i][0][0]->particleMassFlux.fluxR2 = 0;
			}
			double deltaM = massFlux2 + massFlux1 - bins[i][0][0]->particleMassFlux.fluxR1 - bins[i][0][0]->particleMassFlux.fluxR2;
			double deltaP = momentaFlux2 + momentaFlux1 - bins[i][0][0]->particleMomentaFlux.fluxR1 - bins[i][0][0]->particleMomentaFlux.fluxR2;
			//TODO знак U!
			double p = bins[i][0][0]->density*bins[i][0][0]->volume*bins[i][0][0]->U/(sqrt(1 - sqr(bins[i][0][0]->U/speed_of_light)));
			//if(abs(deltaP) > abs(p)){
				//printf("aaaaaaaaa");
			//}
			p = p + deltaP;
			double m = bins[i][0][0]->density*bins[i][0][0]->volume;
			bins[i][0][0]->density += deltaM/(bins[i][0][0]->volume);
			/*if(bins[i][0][0]->density < 0){
				if(abs(bins[i][0][0]->density) < epsilon){
					bins[i][0][0]->density = 0;
				} else {
					printf("aaa");
				}
			}*/
			bins[i][0][0]->massVelocity = p/sqrt(sqr(m) + sqr(p/speed_of_light));
			if(abs(bins[i][0][0]->massVelocity) > speed_of_light){
				if(abs(bins[i][0][0]->massVelocity) < (1 + epsilon)*speed_of_light){
					bins[i][0][0]->massVelocity = (1 - epsilon)*speed_of_light;
				}
				printf("aaa");
			}
		}
	}

	std::list<Particle*>::iterator it = introducedParticles.begin();
	while(it != introducedParticles.end()){
		Particle* particle = *it;
		int* index = SpaceBin::binByCoordinates(particle->absoluteZ, particle->getAbsoluteTheta(), particle->getAbsolutePhi(), upstreamR, deltaR, deltaTheta, deltaPhi);
		if((index[0] > 0)&&(index[0] < rgridNumber)&&(index[1] > 0)&&(index[1] < thetagridNumber)&&(index[2] > 0)&&(index[2] < phigridNumber)){
			SpaceBin* bin = bins[index[0]][index[1]][index[2]];
			particle->setLocalMomentum(bin->U, bin->UTheta, bin->UPhi);
		}
		++it;
	}
}

////обнуляет счётчики зарегистрированных частиц
void Simulation::resetDetectors(){
	zeroBin->resetDetectors();
	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < thetagridNumber; ++j){
			for(int k = 0; k < phigridNumber; ++k){
				SpaceBin* bin = bins[i][j][k];
				bin->resetDetectors();			
		/*if(i >= zeroPoint){
			bin->density = density0*Rtot;
		}*/
			}
		}	
	}
}

////возвращает список частиц, размером number. Частицы распределены по максвеллу
std::list <Particle> Simulation::getParticleGaussDistribution(int number){
	std::list <Particle> l = std::list <Particle>();
	for (int i = 0; i < number; ++i){
		double x = uniRandom() - 0.5;
		double y = uniRandom() - 0.5;
		double z = uniRandom() - 0.5;
		double theta = acos(z/sqrt(x*x + y*y +z*z));
		double phi = atan2(y,x);
		Particle particle = Particle(upstreamR*sin(theta)*cos(phi),upstreamR*sin(theta)*sin(phi),upstreamR*cos(theta),temperature,A,Z, U0, theta, phi);
		particle.weight = 1.0/number;
		l.push_front(particle);
	}
	return l;
}
////возвращает функцию распределения максвелла от импульса, температуры и массы.
double Simulation::maxwell(double p, double temperature, double mass){
	return p*p*exp(-p*p/(2*mass*kBoltzman*temperature))*1/(sqrt(2*pi*mass*kBoltzman*temperature)*sqrt(2*pi*mass*kBoltzman*temperature)*sqrt(2*pi*mass*kBoltzman*temperature));
}
////пересчитывает спектральную плотность магнитного поля в соответствии с уравнением 3.86 у Владимирова
void Simulation::updateMagneticField(){
}

////градиент скорости
vector3d Simulation::gradientSpeed(int i, int j, int k){
	return vector3d(0,0,0);
}
////производная члена уравнения 3.86, отвечающего за каскад по плотности энергии
double Simulation::cascadingDerivativeW(double w,double k,double rho){
	if (kolmogorovCascading){
		return (3.0/2.0)*power(w,1.0/2.0)*power(k,5.0/3.0)*power(rho,-1.0/2.0);
	} else {
		return 0;
	}
}

////производная плотности энергии по волновому числу. i - номер ячейки по x, j - по k
double Simulation::derivativeFieldK(double* field,int j){
	double kstep = (maxK-minK)/(kgridNumber-1);
	if(j == 0){
		return (field[1]-field[0])/kstep;
	} else {
		return (field[j]-field[j-1])/kstep;
	}
}

////Член уравнения 3.86, отвечающий за нестабильность (резонансная или белловская)
double Simulation::instability(double w,int i,int j){
	return 0;
}

////Производная члена уравнения 3.86, отвечающая за каскад по волновому числу
double Simulation::cascadingDerivativeK(double w,double k,double rho){
	if(kolmogorovCascading){
		double result =  (5.0/3.0)*power(w,3.0/2.0)*power(k,2.0/3.0)*power(rho,-1.0/2.0);
		double x = result*0;
		if((x != x)||(result != result)){
			printf("%s","NaN cDK");
		}
		return result;
	} else {
		return 0;
	}
}

////Член уравнения 3.86, отвечающий за диссипацию, определяемый по формуле 3.76
double Simulation::dissipation(double w, int i, double k){
	return 0;
}


void Simulation::updateMaxMinK(){
}

void Simulation::updateMaxMinP(double& minP, double& maxP){
	maxP =bins[0][0][0]->getMaxP();
	Particle* particle = getAnyParticle();
	if(particle == NULL){
		minP = 1;
	} else {
		minP = particle->absoluteMomentum;
	}
	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < thetagridNumber; ++j){
			for(int k = 0; k < phigridNumber; ++k){
				if(maxP < bins[i][j][k]->getMaxP()){
					maxP = bins[i][j][k]->getMaxP();
					if(minP == 0){
						minP = maxP;
					}
				}
				if(minP > bins[i][j][k]->getMinP()){
					minP = bins[i][j][k]->getMinP();
				}
			}
		}
	}
}

void Simulation::evaluateMagneticField(double* startField,double* endField,double deltax,int gridNumber,double gradU, double density,double U,int binNumber){
}

void Simulation::updatePressureSpectralDensity(){
}

std::list <Particle*> Simulation::getParticles(){
	std::list<Particle*> list = std::list<Particle*>();
	int cosmicRayNumber = 0;
	int particleBinNumber = 0;
	allParticlesNumber = 0;
	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < thetagridNumber; ++j){
			for(int k = 0; k < phigridNumber; ++k){
				particleBinNumber = 0;
				SpaceBin* bin = bins[i][j][k];
				bin->initialMomentum = 0;
				//double c2 = speed_of_light*speed_of_light;
				for( int l = 0; l < particlesNumber; ++l){
					//if(order(bin->phi1, phi, bin->phi2) && order(bin->theta1, theta, bin->theta2) && order(bin->r1, r, bin->r2)){
						++allParticlesNumber;
						printf("%d",allParticlesNumber);
						printf("%s","\n");
						Particle* particle = new Particle( A, Z,bin, true);
						particle->weight /= particlesNumber;
						startPDF.push_front(particle);
						list.push_front(particle);
						double v = particle->getAbsoluteV();
						double vr = v*cos(particle->absoluteMomentumTheta); 
						bin->initialMomentum += vr*particle->mass*particle->weight/sqrt(1 - (v*v)/c2);
					//}
				}
			}
		}
	}
	return list;
}

SpaceBin* Simulation::getStartBin(double theta, double phi){
	double deltaTheta = pi/thetagridNumber;
	double deltaPhi = 2*pi/phigridNumber;
	int j = lowerInt(theta / deltaTheta);
	int k = lowerInt(phi / deltaPhi);
	return bins[0][j][k];
}

void Simulation::detectFromTo(int fromIndexR, int fromIndexTheta, int fromIndexPhi, int toIndexR, int toIndexTheta, int toIndexPhi, const Particle& particle){
}

void Simulation::multiplyParticleWeights(double v){
	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < thetagridNumber; ++j){
			for(int k = 0; k < phigridNumber; ++k){
				bins[i][j][k]->multiplyParticleWeights(v);
			}
		}
	}
}

Particle* Simulation::getAnyParticle(){
	for(int j = 0; j < thetagridNumber; ++j){
		for(int k = 0; k < phigridNumber; ++k){
			if(bins[0][j][k]->detectedParticlesR2.size() > 0){
				Particle* particle = *bins[0][j][k]->detectedParticlesR2.begin();
				return particle;
			}
		}
	}
	return NULL;
}

void Simulation::introduceNewParticles(){
	double R = upstreamR - deltaR/2;
	double Theta = deltaTheta/2;
	double Phi = deltaPhi/2;
	std::list<Particle*> list = std::list<Particle*>();

	zeroBin->initialMomentum = 0;
	for( int l = 0; l < particlesNumber; ++l){
		Particle* particle = new Particle( A, Z,zeroBin, true);
		particle->weight /= particlesNumber;
		startPDF.push_front(particle);
		list.push_front(particle);
		double v = particle->getAbsoluteV();
		double vr = v*cos(particle->absoluteMomentumTheta); 
		zeroBin->initialMomentum += vr*particle->mass*particle->weight/sqrt(1 - (v*v)/c2);
	}

	
	std::list<Particle*>::iterator it = list.begin();
	while( it != list.end()){
		Particle* particle = *it;
		++it;
		bool side = false;
		double r = particle->getAbsoluteR();
		double theta = acos(particle->absoluteZ/r);
		double phi =atan2(particle->absoluteY, particle->absoluteX);
		if(phi < 0){
			phi = phi + 2*pi;
		}
		double time = 0;
		int* index = zeroBin->propagateParticle(particle, time, timeStep);
		if(time < timeStep){
			while((index[0] >= 0)){
				SpaceBin* bin = bins[index[0]][index[1]][index[2]];
				if(time != time){
					printf("aaa");
				}
				int* tempIndex = bin->propagateParticle(particle,time, timeStep);
				delete[] index;
				index = tempIndex;
				if(index[0] > rgridNumber){
					printf("aaa");
				}
				if(index[0] == rgridNumber){
					index[0] = rgridNumber - 1;
					bin->detectParticleR2(particle);
					particle->absoluteMomentumTheta = pi - particle->absoluteMomentumTheta;
					particle->absoluteMomentumPhi = 2*pi - particle->absoluteMomentumPhi;
					if(particle->absoluteMomentumPhi > 2*pi){
						particle->absoluteMomentumPhi -= 2*pi;
					}
					particle->moveToBinRight(bin);
				}
				if (time >= timeStep){
					break;
				}
			}
			delete[] index;
		}
	}

	it = zeroBin->detectedParticlesR2.begin();
	while(it != zeroBin->detectedParticlesR2.end()){
		Particle* particle = *it;
		introducedParticles.push_front(new Particle(*particle));
		++it;
	}

	it = list.begin();
	while(it != list.end()){
		Particle* particle = *it;
		delete particle;
		++it;
	}
	list.clear();

}

void Simulation::removeEscapedParticles(){
	std::list<Particle*> list = std::list<Particle*>();
	std::list<Particle*>::iterator it = introducedParticles.begin();
	while(it != introducedParticles.end()){
		Particle* particle = *it;
		++it;
		if(particle->absoluteZ < 0){
			delete particle;
		} else {
			list.push_back(particle);
		}
	}
	introducedParticles.clear();
	introducedParticles = list;
}

void Simulation::collectAverageVelocity(){
	int* count = new int[rgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		count[i] = 0;
		averageVelocity[i] = 0;
	}

	std::list<Particle*>::iterator it = introducedParticles.begin();
	while(it != introducedParticles.end()){
		Particle* particle = *it;
		++it;
		double r = particle->getAbsoluteR();
		double theta = acos(particle->absoluteZ/r);
		double phi =atan2(particle->absoluteY, particle->absoluteX);
		if(phi < 0){
			phi = phi + 2*pi;
		}
		int* index = SpaceBin::binByCoordinates(particle->absoluteZ, theta, phi,upstreamR,deltaR,deltaTheta,deltaPhi);
		if(index[0] >= 0){
			count[index[0]] += 1;
			averageVelocity[index[0]] += particle->getAbsoluteV()*cos(particle->absoluteMomentumTheta);
		}
		delete[] index;
	}

	for(int i = 0; i < rgridNumber; ++i){
		averageVelocity[i] /= count[i];
		bins[i][0][0]->averageVelocity = averageVelocity[i];
		bins[i][0][0]->U = averageVelocity[i];
	}

	delete[] count;
}

void Simulation::resetVelocity(){
	std::list<Particle*>* particles = new std::list<Particle*>[rgridNumber];


	std::list<Particle*>::iterator it = introducedParticles.begin();
	while(it != introducedParticles.end()){
		Particle* particle = *it;
		++it;
		double r = particle->getAbsoluteR();
		double theta = acos(particle->absoluteZ/r);
		double phi =atan2(particle->absoluteY, particle->absoluteX);
		if(phi < 0){
			phi = phi + 2*pi;
		}
		int* index = SpaceBin::binByCoordinates(particle->absoluteZ, theta, phi,upstreamR,deltaR,deltaTheta,deltaPhi);
		if(index[0] >= 0){
			particles[index[0]].push_back(new Particle(*particle));
		}
		delete[] index;
	}

	for(int i = 0; i < rgridNumber; ++i){
		bins[i][0][0]->resetVelocity(particles[i]);
		it = particles[i].begin();
		while( it != particles[i].end()){
			Particle* particle = *it;
			delete particle;
			++it;
		}
		particles[i].clear();
	}

	delete[] particles;
}

int Simulation::maxVelocityDerivativeIndex(){
	int result = 1;
	double maxDerivative = (bins[2][0][0]->U - bins[0][0][0]->U);
	for(int i = 2; i < rgridNumber - 1; ++i){
		//< because derivative must be < 0
		if(bins[i + 1][0][0]->U - bins[i - 1][0][0]->U < maxDerivative) {
			maxDerivative = bins[i+1][0][0] - bins[i - 1][0][0];
			result = i;
		}
	}
	return result;
}

double Simulation::secondVelocityDerivative(int i){
	return (bins[i - 1][0][0]->U + bins[i + 1][0][0]->U - 2*bins[i][0][0]->U)/(deltaR*deltaR);
}

double Simulation::thirdVelocityDerivative(int i){
	return (secondVelocityDerivative(i) - secondVelocityDerivative(i - 1))/deltaR;
}

double Simulation::solveSecondDerivativeZero(int i){
	double derivative2 = secondVelocityDerivative(i);
	double derivative3 = thirdVelocityDerivative(i);
	return (bins[i][0][0]->r - derivative2/derivative3);
}

void Simulation::findShockWavePoint(){
	int i = maxVelocityDerivativeIndex();
	double nextPoint = solveSecondDerivativeZero(i);
	shockWavePoints.push_back(nextPoint);
	shockWaveVelocity.push_back((nextPoint - prevPoint)/timeStep);
	prevPoint = nextPoint;
}