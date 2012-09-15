#include "stdafx.h"
#include "simulation.h"
#include "SpaceBin.h"
#include "particle.h"
#include "output.h"
#include "util.h"
#include <omp.h>
//#include <Windows.h>

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
	energy = 0;
	theorEnergy = 0;
	momentumX = 0;
	theorMomentumX = 0;
	A = 1;
	Z = 1;
	simulationType = 1;
	startPDF = std::list <Particle*>();
	timeStep = defaultTimeStep;
	zeroBin = NULL;
}

Simulation::~Simulation(){
	kolmogorovCascading = true;
	resonantInstability = true;
	bellInstability = true;
	alpha = 1;
	beta = 1;
	gamma = 1;
	delta = 1;
	delete[] averageVelocity;
}

void Simulation::initializeProfile(){
	averageVelocity = new double[rgridNumber];
 	minK = defaultMinK;
	maxK = defaultMaxK;
	deltaR = (downstreamR - upstreamR)/(rgridNumber );
	deltaTheta = pi/(thetagridNumber );
	deltaPhi = 2*pi/(phigridNumber );
	double R = upstreamR + deltaR/2;
	double Theta = deltaTheta/2;
	double Phi = deltaPhi/2;
	bins = new SpaceBin*[rgridNumber];
	zeroBinScale = 2*U0*defaultTimeStep/deltaR;
	zeroBin = new SpaceBin(-zeroBinScale*deltaR/2,zeroBinScale*deltaR,U0,density0,temperature,B0,-1, smallAngleScattering);
    for(int i = 0; i < rgridNumber; ++i){
		double density = density0;
		double u;
		if(i < shockWavePoint){
			u = U0;
		} else {
			u = U0;
		}
		bins[i] = new SpaceBin(R,deltaR,u,density,temperature,B0,i, smallAngleScattering);
		averageVelocity[i] = U0;
		R = R + deltaR;
	}
}
/////   .
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
		//printf("%s", "\n");
		int j = 0;
		int l = 0;
		//printf("%s", "Iteration started\n");
		if(itNumber == 0){
			printf("%s","First iteration\n");
		} else {
			//printf("%s", "Particle propagation\n");
			int it;
			#pragma omp parallel for private(it) shared(l)
			for(it = 0; it < introducedParticles.size(); ++it){
				Particle* particle = introducedParticles[it];
				++l;
				bool side = false;
				double r = particle->absoluteX;

				int index = SpaceBin::binByCoordinates(particle->absoluteX,upstreamR,deltaR, rgridNumber);

				double time = 0;
				while((index >= 0)){
					j++;
					SpaceBin* bin = bins[index];
					if(time != time){
						printf("time != time\n");
					}
					int tempIndex = bin->propagateParticle(particle,time, timeStep, rgridNumber);
					index = tempIndex;
					if(index > rgridNumber){
						printf("index[0] > rgridNumber\n");
					}
					if(index == rgridNumber){
						if((simulationType == 1) || (simulationType == 3)){
							index = rgridNumber - 1;
							particle->absoluteMomentumX= - particle->absoluteMomentumX;
							particle->moveToBinRight(bin);
							theorMomentumX += 2*particle->absoluteMomentumX*particle->weight;
						} else if (simulationType == 2) {
							break;
						}
					}
					if(index <= -1){
						if(simulationType == 3){
							index = 0;
							particle->absoluteMomentumX= - particle->absoluteMomentumX;
							particle->moveToBinLeft(bin);
							theorMomentumX += 2*particle->absoluteMomentumX*particle->weight;
						}
					}
					if (time >= timeStep){
						break;
					}
				}
			}

			//printf("%s", "Introducing new particles from left boundary\n");
			if((simulationType == 1) || (simulationType == 2)){
				introduceNewParticles();
			}
			//printf("%s", "collecting velocity, density and crflux\n");
			collectAverageVelocity();
			//printf("%s", "removing escaped particles\n");
			removeEscapedParticles();
			//printf("%s","updating energy\n");
			updateCosmicRayBoundMomentum(itNumber % 20 == 0);

			outputEnergyPDF(introducedParticles,"./output/tamc_energy_pdf.dat");
			resetDetectors();
			//printf("%s","iteration  ");
			//printf("%d\n",itNumber);
		}
		updateEnergy();
		if(itNumber % 20 == 0){

			if(bins[0]->particles.size() > 0){
				outputPDF(bins[0]->particles,"./output/tamc_pdf0.dat");
			}

			outputEnergyPDF(introducedParticles,"./output/tamc_energy_pdf.dat");
			findShockWavePoint();
			printf("%s","iteration  ");
			printf("%d\n",itNumber);
			printf("outputing\n");
			outputParticles(introducedParticles,"./output/particles.dat");
			outputPDF(introducedParticles,"./output/tamc_pdf.dat");
			outputEnergyPDF(introducedParticles,"./output/tamc_energy_pdf.dat");
			outIteration = fopen("./output/tamc_iteration.dat","a");
			radialFile = fopen("./output/tamc_radial_profile.dat","a");
			fprintf(outIteration,"%d %lf %lf %lf %lf %d %lf %lf\n",itNumber, energy, theorEnergy, momentumX, theorMomentumX, introducedParticles.size(), particlesWeight, shockWavePoint);
			fclose(outIteration);
			outputRadialProfile(bins,radialFile, rgridNumber);
			fclose(radialFile);
		}
	}
}

////   
void Simulation::resetDetectors(){
	zeroBin->resetDetectors();
	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < thetagridNumber; ++j){
			for(int k = 0; k < phigridNumber; ++k){
				SpaceBin* bin = bins[i];
				bin->resetDetectors();			
			}
		}	
	}
}

void Simulation::introduceNewParticles(){
	double R = upstreamR - deltaR/2;
	double Theta = deltaTheta/2;
	double Phi = deltaPhi/2;
	std::vector<Particle*> list = std::vector<Particle*>();

	zeroBin->initialMomentum = 0;
	for( int l = 0; l < zeroBinScale*particlesNumber; ++l){
		Particle* particle = new Particle( A, Z,zeroBin, true, allParticlesNumber);
		allParticlesNumber++;
		particle->weight /= zeroBinScale*particlesNumber;
		list.push_back(particle);
		double v = particle->getAbsoluteV();
		if(abs(v*v/c2 - 1) < epsilon){
			printf("v = c\n in introdeceNewParticles");
		}
		zeroBin->initialMomentum += particle->getAbsoluteVX()*particle->mass*particle->weight/sqrt(1 - (v*v)/c2);
	}

	int i;
	#pragma omp parallel for private(i)
	for(i = 0; i < list.size(); ++i){
		Particle* particle = list[i];
		bool side = false;
		double r = particle->absoluteX;

		double time = 0;
		int index = zeroBin->propagateParticle(particle, time, timeStep, rgridNumber);
		if(time < timeStep){
			while(index >= -1){
				SpaceBin* bin;
				if(index < 0){
					bin = zeroBin;
				} else {
					bin = bins[index];
				}
				if(time != time){
					printf("time != time\n");
				}
				int tempIndex = bin->propagateParticle(particle,time, timeStep, rgridNumber);
				index = tempIndex;
				if(index > rgridNumber){
					printf("index > rgridNumber\n");
				}
				if(index == rgridNumber){
					index = rgridNumber - 1;
					bin->detectParticleR2(particle);
					particle->absoluteMomentumX = - particle->absoluteMomentumX;
					particle->moveToBinRight(bin);
				}
				if (time >= timeStep){
					break;
				}
			}
		}
	}

	std::vector<Particle*>::iterator it = list.begin();
	while(it != list.end()){
		Particle* particle = *it;
		if(particle->absoluteX > 0){
			introducedParticles.push_back(new Particle(*particle));
			theorEnergy += particle->getEnergy()*particle->weight;
			theorMomentumX += particle->absoluteMomentumX*particle->weight;
		}
		delete particle;
		++it;
	}
	list.clear();

}

////  ,  number.    
std::list <Particle> Simulation::getParticleGaussDistribution(int number){
	std::list <Particle> l = std::list <Particle>();
	for (int i = 0; i < number; ++i){
		double x = uniRandom() - 0.5;
		double y = uniRandom() - 0.5;
		double z = uniRandom() - 0.5;
		double theta = acos(z/sqrt(x*x + y*y +z*z));
		if( abs(x*x + y*y +z*z) < epsilon){
			theta = pi/2;
		}
		double phi = atan2(y,x);
		Particle particle = Particle(upstreamR*sin(theta)*cos(phi),upstreamR*sin(theta)*sin(phi),upstreamR*cos(theta),temperature,A,Z, U0, theta, phi);
		particle.weight = 1.0/number;
		l.push_front(particle);
	}
	return l;
}
////     ,   .
double Simulation::maxwell(double p, double temperature, double mass){
	return p*p*exp(-p*p/(2*mass*kBoltzman*temperature))*1/(sqrt(2*pi*mass*kBoltzman*temperature)*sqrt(2*pi*mass*kBoltzman*temperature)*sqrt(2*pi*mass*kBoltzman*temperature));
}
////         3.86  
void Simulation::updateMagneticField(){
}

//// 
vector3d Simulation::gradientSpeed(int i, int j, int k){
	return vector3d(0,0,0);
}

void Simulation::updateMaxMinK(){
}

void Simulation::updateMaxMinP(double& minP, double& maxP){
	maxP =bins[0]->getMaxP();
	Particle* particle = getAnyParticle();
	if(particle == NULL){
		minP = 1;
	} else {
		minP = particle->absoluteMomentum;
	}
	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < thetagridNumber; ++j){
			for(int k = 0; k < phigridNumber; ++k){
				if(maxP < bins[i]->getMaxP()){
					maxP = bins[i]->getMaxP();
					if(minP == 0){
						minP = maxP;
					}
				}
				if(minP > bins[i]->getMinP()){
					minP = bins[i]->getMinP();
				}
			}
		}
	}
}

void Simulation::evaluateMagneticField(double* startField,double* endField,double deltax,int gridNumber,double gradU, double density,double U,int binNumber){
}


std::vector <Particle*> Simulation::getParticles(){
	std::vector<Particle*> list = std::vector<Particle*>();
	int cosmicRayNumber = 0;
	int particleBinNumber = 0;
	allParticlesNumber = 0;
	for(int i = 0; i < rgridNumber; ++i){
		for(int j = 0; j < thetagridNumber; ++j){
			for(int k = 0; k < phigridNumber; ++k){
				particleBinNumber = 0;
				SpaceBin* bin = bins[i];
				bin->initialMomentum = 0;
				int l;
				for(l = 0; l < particlesNumber; ++l){
					++allParticlesNumber;
					printf("%d",allParticlesNumber);
					printf("%s","\n");
					Particle* particle = new Particle( A, Z,bin, true, allParticlesNumber);
					particle->weight /= particlesNumber;
					startPDF.push_front(particle);
					list.push_back(particle);
					double v = particle->getAbsoluteV(); 
					if(abs(v*v/c2 - 1) < epsilon){
						printf("v = c\n in getParticles");
					}
					bin->initialMomentum += particle->absoluteMomentumX*particle->weight;
					energy += particle->getEnergy()*particle->weight;
					momentumX += particle->absoluteMomentumX*particle->weight;

				}
			}
		}
	}
	theorEnergy = energy;
	theorMomentumX = momentumX;
	return list;
}

SpaceBin* Simulation::getStartBin(double theta, double phi){
	double deltaTheta = pi/thetagridNumber;
	double deltaPhi = 2*pi/phigridNumber;
	int j = lowerInt(theta/deltaTheta);
	if( theta == pi ){
		theta = thetagridNumber - 1;
	}
	int k = lowerInt(phi/deltaPhi);
	if (phi == 2*pi){
		k = 0;
	}
	return bins[0];
}

void Simulation::detectFromTo(int fromIndexR, int fromIndexTheta, int fromIndexPhi, int toIndexR, int toIndexTheta, int toIndexPhi, const Particle& particle){
}

Particle* Simulation::getAnyParticle(){
	for(int j = 0; j < thetagridNumber; ++j){
		for(int k = 0; k < phigridNumber; ++k){
			if(bins[0]->detectedParticlesR2.size() > 0){
				Particle* particle = *bins[0]->detectedParticlesR2.begin();
				return particle;
			}
		}
	}
	return NULL;
}

void Simulation::removeEscapedParticles(){
	std::vector<Particle*> list = std::vector<Particle*>();
	std::vector<Particle*>::iterator it = introducedParticles.begin();
	while(it != introducedParticles.end()){
		Particle* particle = *it;
		++it;
		if((particle->absoluteX < 0) || (particle->absoluteX > downstreamR)){
			theorEnergy -= particle->getEnergy()*particle->weight;
			theorMomentumX -= particle->absoluteMomentumX*particle->weight;
			delete particle;
		} else {
			if(particle->absoluteMomentum > momentumParameter*particle->previousAbsoluteMomentum){
				for(int i = 0; i < generationSize; ++i){
					Particle* particle1 = new Particle(*particle);
					particle1->weight /= generationSize;
					particle1->previousAbsoluteMomentum = particle1->absoluteMomentum;
					list.push_back(particle1);
				}
				delete particle;
			} else {
				if(momentumParameter*particle->absoluteMomentum < particle->previousAbsoluteMomentum){
					particle->previousAbsoluteMomentum = particle->absoluteMomentum;
				}
				list.push_back(particle);
			}
		}
	}
	introducedParticles.clear();
	introducedParticles = list;
}

void Simulation::collectAverageVelocity(){
	double* count = new double[rgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		count[i] = 0;
		averageVelocity[i] = 0;
		bins[i]->density = 0;
	}

	std::vector<Particle*>::iterator it = introducedParticles.begin();

	while(it != introducedParticles.end()){
		Particle* particle = *it;
		++it;
		double theta;
		double r = particle->absoluteX;

		int index = SpaceBin::binByCoordinates(particle->absoluteX,upstreamR,deltaR, rgridNumber);
		if((index >= 0)&&(index < rgridNumber)){
			count[index] += particle->weight;
			averageVelocity[index] += particle->getAbsoluteVX()*particle->weight;
			bins[index]->density += particle->mass*particle->weight;
			bins[index]->particleMomentaZ.push_back(particle->localMomentumX);
			bins[index]->particleWeights.push_back(particle->weight);
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		bins[i]->density /= bins[i]->volume;
		if(count[i] > epsilon){
			averageVelocity[i] /= count[i];
			bins[i]->averageVelocity = averageVelocity[i];
			bins[i]->U = averageVelocity[i];
		} else {
			printf("0 particles in bin\n");
			bins[i]->averageVelocity = 0;
			bins[i]->U = 0;
		}
	}

	delete[] count;
}

void Simulation::sortParticlesIntoBins(){
	std::vector<Particle*>::iterator it = introducedParticles.begin();\
	for(int i = 0; i < rgridNumber; ++i){
		bins[i]->density = 0.0;
	}
	while( it != introducedParticles.end()){
		Particle* particle = *it;
		double r = particle->absoluteX;

		int index = SpaceBin::binByCoordinates(particle->absoluteX, upstreamR, deltaR, rgridNumber);
		if((index >= 0) && (index < rgridNumber)){
			//bins[index[0]][index[1]][index[2]]->particles.push_back(new Particle(*particle));
			bins[index]->particles.push_back(particle);
			bins[index]->density += particle->mass*particle->weight;
		}
		++it;
	}
	for(int i = 0; i < rgridNumber; ++i){
		bins[i]->density /= bins[i]->volume;
	}
}

void Simulation::updateCosmicRayBoundMomentum(bool write){
	for(int i = 0; i < rgridNumber; ++i){
		bins[i]->updateCosmicRayBoundMomentum(write);
	}
}

void Simulation::smoothProfile(){
	printf("smoothing profile\n");
	for(int i = 0; i < rgridNumber - 1; ++i){
		//if((abs(bins[i][0][0]->U - bins[i + 1][0][0]->U) > 0.4*abs(bins[i][0][0]->U + bins[i + 1][0][0]->U)) || (abs(bins[i][0][0]->density - bins[i + 1][0][0]->density) > 0.4*abs(bins[i][0][0]->density + bins[i + 1][0][0]->density))){
		if(abs(bins[i]->U - bins[i + 1]->U) > 0.3*abs(bins[i]->U + bins[i + 1]->U)){
			std::list<SpaceBin*> list;
			list.push_back(bins[i]);
			list.push_back(bins[i + 1]);
			smoothProfile(list);
			list.clear();
		}
	}
}

void Simulation::smoothProfile(std::list<SpaceBin*> bins){
	std::list<SpaceBin*>::iterator it = bins.begin();

	double volume = 0;
	double u = 0;
	double mass = 0;
	double count = 0;
	double density = 0;

	while(it != bins.end()){
		SpaceBin* bin = *it;
		volume += bin->volume;
		
		mass += bin->density*bin->volume;

		count += bin->particles.size();

		u += bin->U*bin->particles.size();

		++it;
	}


	if( count > 0){
		u /= count;
	} else {
		u = 0;
	}

	if(bins.size() > 0){
		density = mass/volume;
	} else {
		density = 0;
	}

	it = bins.begin();
	while(it != bins.end()){
		SpaceBin* bin = *it;
		
		bin->U = u;
		bin->density = density;

		++it;
	}
}

void Simulation::updateEnergy(){
	std::vector<Particle*>::iterator it = introducedParticles.begin();
	energy = 0;
	momentumX = 0;
	particlesWeight = 0;
	while(it != introducedParticles.end()){
		Particle* particle = *it;
		energy += particle->getEnergy()*particle->weight;
		if(particle->absoluteMomentum < 0){
			printf("particle->absoluteMomentum < 0\n");
		}
		momentumX += particle->absoluteMomentumX*particle->weight;
		particlesWeight += particle->weight;
		++it;
	}
}

void Simulation::findShockWavePoint(){
	double deltaV = 0;
	int shockWaveIndex = 0;
	for(int i = 1; i < rgridNumber - 1; ++i){
		double dV = bins[i - 1]->U - bins[i]->U;
		if( dV > deltaV){
			shockWaveIndex = i;
			deltaV = dV;
		}
	}
	shockWavePoint = bins[shockWaveIndex]->r;
}