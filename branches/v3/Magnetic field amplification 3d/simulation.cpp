#include "stdafx.h"
#include "simulation.h"
#include "SpaceBin.h"
#include "particle.h"
#include "output.h"
#include "util.h"
#include "constants.h"
#include <omp.h>
//#include <Windows.h>

Simulation::Simulation(){
	partOfCosmicRay = 0;
	allParticlesNumber = 0;
	kolmogorovCascading = true;
	resonantInstability = true;
	bellInstability = true;
	epsilonR = 0.1;
	A = 1;
	Z = 1;
	startPDF = std::list <Particle*>();
	timeStep = defaultTimeStep;
	shockWavePoint = 0;
	theorEnergy = 0;
	theorMomentumZ = 0;
	theorMomentumY = 0;
	theorMomentumX = 0;
	energy = 0;
	momentumZ = 0;
	momentumY = 0;
	momentumX = 0;
}

Simulation::~Simulation(){
	kolmogorovCascading = true;
	resonantInstability = true;
	bellInstability = true;
	epsilonR = 0.1;
	delete[] averageVelocity;
	/*for(int i = 0; i < xgridNumber; ++i){
		delete[] pressureSpectralDensity[i];
	}
	delete[] pressureSpectralDensity;*/
}

void Simulation::initializeProfile(){
	averageVelocity = new double[rgridNumber];
	currentShockWavePoint = shockWavePoint;
 	minK = defaultMinK;
	maxK = defaultMaxK;
	deltaR = (downstreamR - upstreamR)/(rgridNumber );
	deltaTheta = pi/(thetagridNumber );
	deltaPhi = 2*pi/(phigridNumber );
	double R = upstreamR + deltaR/2;
	double Theta = deltaTheta/2;
	double Phi = deltaPhi/2;
	bins = new SpaceBin***[rgridNumber];
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
					//u = U0*sqr((upstreamR + deltaR/2)/R);
					u = U0;
					//density = density0*sqr((upstreamR + deltaR/2)/R);
				} else {
					u = U0/R;
					//u = U0*sqr((upstreamR + deltaR/2)/R)/Rtot;
					//density = density0*sqr((upstreamR + deltaR/2)/R)/Rtot;
				}
				//u = U0;
				bins[i][j][k] = new SpaceBin(R,Theta,Phi,deltaR,deltaTheta,deltaPhi,u,density,Theta,Phi,temperature,B0,i,j,k,smallAngleScattering);
				Phi = Phi + deltaPhi;
			}
			Theta = Theta + deltaTheta;
		}
		R = R + deltaR;
	}
	/*for(int j = 0; j < xgridNumber; ++j){
		for(int i =0; i < kgridNumber; ++i){
			xbins[j]->minK = defaultMinK;
			xbins[j]->maxK = defaultMaxK;
			double k = minK + (maxK-minK)*i/(kgridNumber-1);
			xbins[j]->magneticField[i] = turbulenceSeed*turbulenceSeed/(4*pi*maxK*k*log(maxK/minK));
			xbins[j]->updateMagneticField();
		}
	}*/
}
///// ������� ������� ���������.
void Simulation::simulate(){
	FILE* outIteration;
	srand ( time(NULL) );
	initializeProfile();
	//FILE* outIteration = fopen("./output/tamc_iteration.dat","w");
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
			printf("%s","First iteration\n");
		} else {
			printf("%s", "Particle propagation\n");
			int it;
			if(introducedParticles.size() == 0){
				break;
			}
			#pragma omp parallel for private(it) shared(l)
			for(it = 0; it < introducedParticles.size(); ++it){
				Particle* particle = introducedParticles[it];
				//printf("%d %s",l,"\n");
				++l;
				bool side = false;
				double r = sqrt(particle->absoluteX*particle->absoluteX + particle->absoluteY*particle->absoluteY + particle->absoluteZ*particle->absoluteZ);
				double theta = acos(particle->absoluteZ/r);
				if(abs(r) < epsilon){
					theta = pi/2;
				}
				double phi =atan2(particle->absoluteY, particle->absoluteX);
				if(phi < 0){
					phi = phi + 2*pi;
				}
				int* index = SpaceBin::binByCoordinates(r, theta, phi,upstreamR,deltaR,deltaTheta,deltaPhi, rgridNumber);
				double time = 0;
				while((index[0] >= 0)&&(index[0] < rgridNumber)){
					j++;
					//if(j > 100){
						//printf("%d \n",number);
					//}
					//printf("%d",number);
					//printf("%s","\n");
					//cout>>number>>"\n";
					//std::cout<<number<<"\n";
					SpaceBin* bin = bins[index[0]][index[1]][index[2]];
					if(time != time){
						printf("time != time\n");
					}
					int* tempIndex = bin->propagateParticle(particle,time, timeStep, rgridNumber);
					delete[] index;
					index = tempIndex;
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
			//printf("%s", "Reseting profile\n");
			//resetProfile();
			printf("%s", "Removing escaped particles\n");
			removeEscapedParticles();
			//Sleep(5000);
			printf("%s", "collect average velocity\n");
			collectAverageVelocity();
			//Sleep(5000);
			//printf("%s", "sorting particles into bins\n");
			//sortParticlesIntoBins();
			//Sleep(5000);
			//printf("%s", "magnetic Field updating\n");
			//updateMagneticField();
			//output(*this);
			updateCosmicRayBoundMomentum();
			printf("%s", "reseting detectors\n");
			resetDetectors();
			//Sleep(5000);
			printf("%s","iteration � ");
			printf("%d\n",itNumber);
		}
		if(introducedParticles.size() > 0){
			updateEnergy();
			updateShockWavePoint();
			if(itNumber % 10 == 0){
				printf("%s", "outputing\n");
				outputParticles(introducedParticles,"./output/particles.dat");
				outputPDF(introducedParticles,"./output/tamc_pdf.dat");
				outputEnergyPDF(introducedParticles,"./output/tamc_energy_pdf.dat");
				outIteration = fopen("./output/tamc_iteration.dat","a");
				radialFile = fopen("./output/tamc_radial_profile.dat","a");
				fprintf(outIteration,"%d %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf\n",itNumber, energy, theorEnergy, momentumZ, theorMomentumZ, momentumX, theorMomentumX, momentumY, theorMomentumY, introducedParticles.size(), particlesWeight, bins[currentShockWavePoint][0][0]->r);
				fclose(outIteration);
				outputRadialProfile(bins,0,0,radialFile, rgridNumber);
				//outputShockWave(shockWavePoints, shockWaveVelocity);
				fclose(radialFile);
			}
			//Sleep(5000);
		}
	}
}

void Simulation::resetProfile(){
	/*for(int i = 0; i < rgridNumber; ++i){
		for( int j = 0; j < thetagridNumber; ++j){
			for(int k = 0; k < phigridNumber; ++k){
				bins[i][j][k]->updateFluxes();
			}
		}
	}*/
	if((thetagridNumber == 1) && (phigridNumber == 1)){
		//TODO what do with i = 0
		for(int i = 0; i < rgridNumber; i++){
			double massFlux1 = 0;
			double momentaFlux1 = 0;
			if(i > 0){
				massFlux1 = bins[i - 1][0][0]->particleMassFlux.fluxR2;
				momentaFlux1 = bins[i - 1][0][0]->particleMomentaFlux.fluxR2;
			}
			double massFlux2 = 0;
			double momentaFlux2 = 0;
			if ( i < rgridNumber - 1){
				massFlux2 = bins[i + 1][0][0]->particleMassFlux.fluxR1;
				momentaFlux2 = bins[i + 1][0][0]->particleMomentaFlux.fluxR1;
			}
			double deltaM = massFlux2 + massFlux1 - bins[i][0][0]->particleMassFlux.fluxR1 - bins[i][0][0]->particleMassFlux.fluxR2;
			double deltaP = momentaFlux2 + momentaFlux1 - bins[i][0][0]->particleMomentaFlux.fluxR1 - bins[i][0][0]->particleMomentaFlux.fluxR2;
			//TODO ���� U!
			/*double p = bins[i][0][0]->density*bins[i][0][0]->volume*bins[i][0][0]->U/(sqrt(1 - sqr(bins[i][0][0]->U/speed_of_light)));
			p = p + deltaP;
			if(abs(deltaP) > abs(p)){
				printf("deltaP > p\n");
			}
			double m = bins[i][0][0]->density*bins[i][0][0]->volume;*/
			bins[i][0][0]->density += deltaM/(bins[i][0][0]->volume);
		}
	}
}

////�������� �������� ������������������ ������
void Simulation::resetDetectors(){
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

////���������� ������ ������, �������� number. ������� ������������ �� ���������
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
////���������� ������� ������������� ��������� �� ��������, ����������� � �����.
double Simulation::maxwell(double p, double temperature, double mass){
	return p*p*exp(-p*p/(2*mass*kBoltzman*temperature))*1/(sqrt(2*pi*mass*kBoltzman*temperature)*sqrt(2*pi*mass*kBoltzman*temperature)*sqrt(2*pi*mass*kBoltzman*temperature));
}
////������������� ������������ ��������� ���������� ���� � ������������ � ���������� 3.86 � �����������

////�������� ��������

////����������� ����� ��������� 3.86, ����������� �� ������ �� ��������� �������
double Simulation::cascadingDerivativeW(double w,double k,double rho){
	if (kolmogorovCascading){
		return (3.0/2.0)*power(w,1.0/2.0)*power(k,5.0/3.0)*power(rho,-1.0/2.0);
	} else {
		return 0;
	}
}

////����������� ��������� ������� �� ��������� �����. i - ����� ������ �� x, j - �� k
double Simulation::derivativeFieldK(double* field,int j){
	double kstep = (maxK-minK)/(kgridNumber-1);
	if(j == 0){
		return (field[1]-field[0])/kstep;
	} else {
		return (field[j]-field[j-1])/kstep;
	}
}

////���� ��������� 3.86, ���������� �� �������������� (����������� ��� ����������)

////����������� ����� ��������� 3.86, ���������� �� ������ �� ��������� �����

////���� ��������� 3.86, ���������� �� ����������, ������������ �� ������� 3.76



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
	/*std::list <Particle> ::iterator it = xbins[0]->detectedParticlesRight.begin();
	while(it != xbins[0]->detectedParticlesRight.end()){
		Particle particle = *it;
		if(particle.momentum < minP){
			minP = particle.momentum;
		}
		if(particle.momentum > maxP){
			maxP = particle.momentum;
		}
		++it;
	}
	for(int i = 0; i< xgridNumber; ++i){
		double p1 = xbins[i]->getMaxP();
		double p2 = xbins[i]->getMinP();
		if(p1 > maxP){
			maxP = p1;
		}
		if(p2 < minP){
			minP = p2;
		}

	}*/
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
				SpaceBin* bin = bins[i][j][k];
				bin->initialMomentum = 0;
				//double c2 = speed_of_light*speed_of_light;
				for( int l = 0; l < particlesNumber; ++l){
					//if(order(bin->phi1, phi, bin->phi2) && order(bin->theta1, theta, bin->theta2) && order(bin->r1, r, bin->r2)){
						++allParticlesNumber;
						printf("%d",allParticlesNumber);
						printf("%s","\n");
						Particle* particle = new Particle( A, Z,bin, true, allParticlesNumber);
						particle->weight /= particlesNumber;
						startPDF.push_front(particle);
						list.push_back(particle);
						double v = particle->getAbsoluteV();
						double vr = particle->getRadialSpeed(); 
						bin->initialMomentum += vr*particle->mass*particle->weight/sqrt(1 - (v*v)/c2);
						theorEnergy += particle->getEnergy()*particle->weight;
						theorMomentumZ += particle->absoluteMomentum*cos(particle->absoluteMomentumTheta)*particle->weight;
						theorMomentumX += particle->absoluteMomentum*sin(particle->absoluteMomentumTheta)*cos(particle->absoluteMomentumPhi)*particle->weight;
						theorMomentumY += particle->absoluteMomentum*sin(particle->absoluteMomentumTheta)*sin(particle->absoluteMomentumPhi)*particle->weight;
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
	int j = lowerInt(theta/deltaTheta);
	if( theta == pi ){
		theta = thetagridNumber - 1;
	}
	int k = lowerInt(phi/deltaPhi);
	if (phi == 2*pi){
		k = 0;
	}
	return bins[0][j][k];
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

void Simulation::collectAverageVelocity(){
	double* weights = new double[rgridNumber];
	int* count = new int[rgridNumber];
	for(int i = 0; i < rgridNumber; ++i){
		weights[i] = 0;
		count[i] = 0;
		averageVelocity[i] = 0;
		bins[i][0][0]->density = 0;
	}

	std::vector<Particle*>::iterator it = introducedParticles.begin();
	while(it != introducedParticles.end()){
		Particle* particle = *it;
		++it;
		double r = particle->getAbsoluteR();
		double theta = acos(particle->absoluteZ/r);
		if( abs(r) < epsilon){
			theta = pi/2;
		}
		double phi =atan2(particle->absoluteY, particle->absoluteX);
		if(phi < 0){
			phi = phi + 2*pi;
		}
		int* index = SpaceBin::binByCoordinates(r, theta, phi,upstreamR,deltaR,deltaTheta,deltaPhi, rgridNumber);
		if(index[0] >= 0){
			weights[index[0]] += particle->weight;
			count[index[0]] += 1;
			averageVelocity[index[0]] += particle->getRadialSpeed()*particle->weight;
			bins[index[0]][index[1]][index[2]]->density += particle->mass*particle->weight;
			bins[index[0]][index[1]][index[2]]->particleMomentaZ.push_back(particle->localMomentumZ);
			bins[index[0]][index[1]][index[2]]->particleWeights.push_back(particle->weight);
		}
		delete[] index;
	}

	for(int i = 0; i < rgridNumber; ++i){
		bins[i][0][0]->density /= bins[i][0][0]->volume;
		if(weights[i] > epsilon){
			averageVelocity[i] /= weights[i];
			if(count[i] < sqrt(1.0*particlesNumber)){
				averageVelocity[i] = 0;
			}
			//bins[i][0][0]->averageVelocity = averageVelocity[i];
			bins[i][0][0]->U = averageVelocity[i];
			bins[i][0][0]->averageVelocity = averageVelocity[i];

		} else {
			bins[i][0][0]->U = 0;
			bins[i][0][0]->averageVelocity = 0;
			printf("0 particles in bin\n");
		}
	}

	delete[] weights;
	delete[] count;
}

void Simulation::sortParticlesIntoBins(){
	for(int i = 0; i < rgridNumber; ++i){
		bins[i][0][0]->density = 0.0;
	}
	std::vector<Particle*>::iterator it = introducedParticles.begin();
	while( it != introducedParticles.end()){
		Particle* particle = *it;
		double r = particle->getAbsoluteR();
		double theta = acos(particle->absoluteZ/r);
		if( abs(r) < epsilon){
			printf("r < epsilon");
			theta = pi/2;
		}
		if( r == 0.0){
			theta = pi/2;
		}
		double phi =atan2(particle->absoluteY, particle->absoluteX);
		if(phi < 0){
			phi = phi + 2*pi;
		}
		int* index = SpaceBin::binByCoordinates(r, theta, phi,upstreamR,deltaR,deltaTheta,deltaPhi, rgridNumber);
		if((index[0] >= 0) && (index[0] < rgridNumber) && (index[1] >= 0) && (index[1] < thetagridNumber) && (index[2] >= 0) && (index[2] < phigridNumber)){
			bins[index[0]][index[1]][index[2]]->particles.push_back(new Particle(*particle));
			bins[index[0]][index[1]][index[2]]->density += particle->mass*particle->weight;
		}
		++it;
	}
	for(int i = 0; i < rgridNumber; ++i){
		bins[i][0][0]->density /= bins[i][0][0]->volume;
	}
}

void Simulation::removeEscapedParticles(){
	std::vector<Particle*> list = std::vector<Particle*>();
	std::vector<Particle*>::iterator it = introducedParticles.begin();
	while(it != introducedParticles.end()){
		Particle* particle = *it;
		++it;
		if(particle->getAbsoluteR() > downstreamR){
			theorEnergy -= particle->getEnergy()*particle->weight;
			theorMomentumZ -= particle->absoluteMomentum*cos(particle->absoluteMomentumTheta)*particle->weight;
			theorMomentumX -= particle->absoluteMomentum*sin(particle->absoluteMomentumTheta)*cos(particle->absoluteMomentumPhi)*particle->weight;
			theorMomentumY -= particle->absoluteMomentum*sin(particle->absoluteMomentumTheta)*sin(particle->absoluteMomentumPhi)*particle->weight;
			delete particle;
		} else {
			if(particle->absoluteMomentum > momentumParameter*particle->previousAbsoluteMomentum){
				for(int i = 0; i < generationSize; ++i){
					Particle* particle1 = new Particle(*particle);
					particle1->weight /= generationSize;
					particle1->previousAbsoluteMomentum = particle1->absoluteMomentum;
					allParticlesNumber++;
					particle1->number = allParticlesNumber;
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

void Simulation::updateEnergy(){
	std::vector<Particle*>::iterator it = introducedParticles.begin();
	energy = 0;
	momentumZ = 0;
	momentumX = 0;
	momentumY = 0;
	particlesWeight = 0;
	while(it != introducedParticles.end()){
		Particle* particle = *it;
		energy += particle->getEnergy()*particle->weight;
		if(particle->absoluteMomentum < 0){
			printf("particle->absoluteMomentum < 0\n");
		}
		momentumZ += particle->absoluteMomentum*cos(particle->absoluteMomentumTheta)*particle->weight;
		momentumX += particle->absoluteMomentum*sin(particle->absoluteMomentumTheta)*cos(particle->absoluteMomentumPhi)*particle->weight;
		momentumY += particle->absoluteMomentum*sin(particle->absoluteMomentumTheta)*sin(particle->absoluteMomentumPhi)*particle->weight;
		particlesWeight += particle->weight;
		++it;
	}
}

void Simulation::updateCosmicRayBoundMomentum(){
	for(int i = 0; i < rgridNumber; ++i){
		bins[i][0][0]->updateCosmicRayBoundMomentum();
	}
}

void Simulation::updateShockWavePoint(){
	double maxDensity = bins[rgridNumber - 1][0][0]->density;
	for(int i = rgridNumber-1;i >=0; --i){
		if(bins[i][0][0]->density > maxDensity){
			maxDensity = bins[i][0][0]->density;
			currentShockWavePoint = i;
		}
	}
	//���������� ������ ������� �����
	double averageDensity = 0;
	double fullVolume = 0;
	for(int i = 0; i < rgridNumber; ++i){
		fullVolume += bins[i][0][0]->volume;
		averageDensity +=bins[i][0][0]->density*bins[i][0][0]->volume;
	}
	averageDensity /= fullVolume;
	if(maxDensity > 2*averageDensity){
		int secondShockWavePoint = currentShockWavePoint;
		int tempShockWavePoint = currentShockWavePoint;

		for(int i = currentShockWavePoint + 1; i < rgridNumber; ++i){
			if((bins[i][0][0]->density - averageDensity) < (maxDensity-averageDensity)/2){
				tempShockWavePoint = i;
				maxDensity = bins[i][0][0]->density;
				break;
			}
		}

		for(int i = tempShockWavePoint; i < rgridNumber; ++i){
			if(bins[i][0][0]->density > maxDensity){
				maxDensity = bins[i][0][0]->density;
				secondShockWavePoint = i;
			}
		}

		if(tempShockWavePoint != secondShockWavePoint){
			currentShockWavePoint = secondShockWavePoint;
		}
	}
}