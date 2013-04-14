#include "stdafx.h"
#include "simulation.h"
#include "SpaceBin.h"
#include "particle.h"
#include "output.h"
#include "util.h"
#include <omp.h>
#include <Windows.h>

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
	deltaT = defaultTimeStep;
	simulationTime = 0;
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
	bins = new SpaceBin*[rgridNumber];
    for(int i = 0; i < rgridNumber; ++i){
		//bins[i] = new SpaceBin**[thetagridNumber];
		Theta = deltaTheta/2;
		//for(int j = 0; j < thetagridNumber; ++j){
			//bins[i][j] = new SpaceBin*[phigridNumber];
			Phi = deltaPhi/2;
			//for(int k = 0; k < phigridNumber; ++k){
			    double density = density0;
				double temperature = temperature0;
				double u;
			    //init shock wave
				if((1 <= i) && (i < 2*shockWavePoint)){
					u = U0;
				} else {
					//u = U0/R;
					u = 0;
				}
				//

				//init sound wave
				  /*double c = sqrt(gamma*gasConstant*temperature0); 
				  density = density0 + 0.01*density0*sin(pi*10.0*i/rgridNumber);
				  u = c*0.01*density0*sin(pi*10.0*i/rgridNumber)/density0;
				  temperature = temperature0 + c*u*(gamma - 1)/(gamma*gasConstant);*/
				//


				bins[i] = new SpaceBin(R,Theta,Phi,deltaR,deltaTheta,deltaPhi,u,density,Theta,Phi,temperature,B0,i,0,0,smallAngleScattering);
				//Phi = Phi + deltaPhi;
			//}
			//Theta = Theta + deltaTheta;
		//}
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
///// главная функция программы.
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
				/*Particle* particle = introducedParticles[it];
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
				delete[] index;*/
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
			printf("%s","iteration № ");
			printf("%d\n",itNumber);
			simulationTime = simulationTime + deltaT;
			printf("%s","time ");
			printf("%lf\n",simulationTime);
		}
		if(introducedParticles.size() > 0){
			updateEnergy();
			updateShockWavePoint();
			if(itNumber % 100 == 0){
				printf("%s", "outputing\n");
				outputParticles(introducedParticles,"./output/particles.dat");
				outputPDF(introducedParticles,"./output/tamc_pdf.dat");
				outputEnergyPDF(introducedParticles,"./output/tamc_energy_pdf.dat");
				outIteration = fopen("./output/tamc_iteration.dat","a");
				radialFile = fopen("./output/tamc_radial_profile.dat","a");
				fprintf(outIteration,"%d %lf %lf %lf %lf %d %lf %lf\n",itNumber, energy, theorEnergy, momentumZ, theorMomentumZ, introducedParticles.size(), particlesWeight, bins[currentShockWavePoint]->r);
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
				massFlux1 = bins[i - 1]->particleMassFlux.fluxR2;
				momentaFlux1 = bins[i - 1]->particleMomentaFlux.fluxR2;
			}
			double massFlux2 = 0;
			double momentaFlux2 = 0;
			if ( i < rgridNumber - 1){
				massFlux2 = bins[i + 1]->particleMassFlux.fluxR1;
				momentaFlux2 = bins[i + 1]->particleMomentaFlux.fluxR1;
			}
			double deltaM = massFlux2 + massFlux1 - bins[i]->particleMassFlux.fluxR1 - bins[i]->particleMassFlux.fluxR2;
			double deltaP = momentaFlux2 + momentaFlux1 - bins[i]->particleMomentaFlux.fluxR1 - bins[i]->particleMomentaFlux.fluxR2;
			//TODO знак U!
			/*double p = bins[i][0][0]->density*bins[i][0][0]->volume*bins[i][0][0]->U/(sqrt(1 - sqr(bins[i][0][0]->U/speed_of_light)));
			p = p + deltaP;
			if(abs(deltaP) > abs(p)){
				printf("deltaP > p\n");
			}
			double m = bins[i][0][0]->density*bins[i][0][0]->volume;*/
			bins[i]->density += deltaM/(bins[i]->volume);
		}
	}
}

////обнуляет счётчики зарегистрированных частиц
void Simulation::resetDetectors(){
	for(int i = 0; i < rgridNumber; ++i){
				SpaceBin* bin = bins[i];
				bin->resetDetectors();			
		/*if(i >= zeroPoint){
			bin->density = density0*Rtot;
		}*/
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
		if( abs(x*x + y*y +z*z) < epsilon){
			theta = pi/2;
		}
		double phi = atan2(y,x);
		//todo!
		Particle particle = Particle(upstreamR*sin(theta)*cos(phi),upstreamR*sin(theta)*sin(phi),upstreamR*cos(theta),temperature0,A,Z, U0, theta, phi);
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

////градиент скорости

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

////Производная члена уравнения 3.86, отвечающая за каскад по волновому числу

////Член уравнения 3.86, отвечающий за диссипацию, определяемый по формуле 3.76



void Simulation::updateMaxMinP(double& minP, double& maxP){
	maxP =bins[0]->getMaxP();
	Particle* particle = getAnyParticle();
	if(particle == NULL){
		minP = 1;
	} else {
		minP = particle->absoluteMomentum;
	}
	for(int i = 0; i < rgridNumber; ++i){
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
				particleBinNumber = 0;
				SpaceBin* bin = bins[i];
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

void Simulation::multiplyParticleWeights(double v){
	for(int i = 0; i < rgridNumber; ++i){
				bins[i]->multiplyParticleWeights(v);
	}
}

Particle* Simulation::getAnyParticle(){
	if(bins[0]->detectedParticlesR2.size() > 0){
				Particle* particle = *bins[0]->detectedParticlesR2.begin();
				return particle;
	}
	return NULL;

}

void Simulation::collectAverageVelocity(){
	/*double* weights = new double[rgridNumber];
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
	delete[] count;*/
	double* newMomentum = new double[rgridNumber];
	double* newDensity = new double[rgridNumber];
	double* newPressure = new double[rgridNumber];

	evaluateHydrodynamic(newDensity, newMomentum, newPressure);

	for(int i = 0; i < rgridNumber; ++i){
		if((newMomentum[i] != newMomentum[i]) || (0*newMomentum[i] != 0*newMomentum[i])){
			printf("NaN velocity\n");
			Sleep(500);
		}

		//bins[i]->U = (newMomentum[i]/bins[i]->density);
		bins[i]->U = (newMomentum[i]/newDensity[i]);

		if(newDensity[i] < 0){
			bins[i]->density = epsilon;
			printf("density < 0\n");
		} else {
			if((newDensity[i] != newDensity[i]) || (0*newDensity[i] != 0*newDensity[i])){
				printf("NaN density\n");
				Sleep(1000);
			}
			bins[i]->density = newDensity[i];
		}

		//bins[i]->U = newMomentum[i]/bins[i]->density;

		if(bins[i]->density < 100*epsilon*density0){
			bins[i]->U = 0;
		}

		//bins[i]->U = newMomentum[i]/bins[i]->density;


		//bins[i]->pressure = (gamma - 1)*(newPressure[i] - bins[i]->density*bins[i]->U*bins[i]->U/2);
		bins[i]->pressure = newPressure[i];

		if(bins[i]->pressure < 0){
			bins[i]->pressure = epsilon;
			printf("pressure < 0\n");
		} else {
			if((bins[i]->pressure != bins[i]->pressure) || (0*bins[i]->pressure != 0*bins[i]->pressure)){
				printf("NaN pressure\n");
				Sleep(500);
			}
		}

		if(bins[i]->density < 100*epsilon*density0){
			bins[i]->pressure = epsilon*newPressure[i];
		}
	}

	delete[] newMomentum;
	delete[] newDensity;
	delete[] newPressure;
}

void Simulation::sortParticlesIntoBins(){
	for(int i = 0; i < rgridNumber; ++i){
		bins[i]->density = 0.0;
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
		int index = SpaceBin::binByCoordinates(r, theta, phi,upstreamR,deltaR,deltaTheta,deltaPhi, rgridNumber);
		if((index >= 0) && (index < rgridNumber)){
			bins[index]->particles.push_back(new Particle(*particle));
			//bins[index]->density += particle->mass*particle->weight;
		}
		++it;
	}
	/*for(int i = 0; i < rgridNumber; ++i){
		bins[i]->density /= bins[i]->volume;
	}*/
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
		bins[i]->updateCosmicRayBoundMomentum();
	}
}

void Simulation::updateShockWavePoint(){
	double maxDensity = bins[rgridNumber - 1]->density;
	for(int i = rgridNumber-1;i >=0; --i){
		if(bins[i]->density > maxDensity){
			maxDensity = bins[i]->density;
			currentShockWavePoint = i;
		}
	}
	//нахождение второй ударной волны
	/*double averageDensity = 0;
	double fullVolume = 0;
	for(int i = 0; i < rgridNumber; ++i){
		fullVolume += bins[i]->volume;
		averageDensity +=bins[i]->density*bins[i]->volume;
	}
	averageDensity /= fullVolume;
	if(maxDensity > 2*averageDensity){
		int secondShockWavePoint = currentShockWavePoint;
		int tempShockWavePoint = currentShockWavePoint;

		for(int i = currentShockWavePoint + 1; i < rgridNumber; ++i){
			if((bins[i]->density - averageDensity) < (maxDensity-averageDensity)/2){
				tempShockWavePoint = i;
				maxDensity = bins[i]->density;
				break;
			}
		}

		for(int i = tempShockWavePoint; i < rgridNumber; ++i){
			if(bins[i]->density > maxDensity){
				maxDensity = bins[i]->density;
				secondShockWavePoint = i;
			}
		}

		if(tempShockWavePoint != secondShockWavePoint){
			currentShockWavePoint = secondShockWavePoint;
		}
	}*/
}

void Simulation::evaluateHydrodynamic(double* newDensity, double* newMomentum, double* newPressure){
	double* tempMomentum = new double[rgridNumber];
	double* tempDensity = new double[rgridNumber];
	double* tempPressure = new double[rgridNumber];
	double c = 2*findMaxVelocity();
	if(abs(c) < epsilon){
		c = scaleParameter*deltaR/defaultTimeStep;
	} 
	deltaT = scaleParameter*deltaR/abs(c);
	//deltaT = min(deltaT, defaultTimeStep);
	//deltaT = min(deltaT, kurzrockMinDeltaT());
	//deltaT = max(deltaT, defaultTimeStep);

	for(int i = 0; i < rgridNumber; ++i){
		//in new U = U + c
		bins[i]->U = bins[i]->U + c;
		newDensity[i] = bins[i]->density;
		newMomentum[i] = bins[i]->density*bins[i]->U;
		newPressure[i] = bins[i]->pressure;
		//newPressure[i] = bins[i]->fullEnergy();
	}

	//LaxVendorf(newDensity, newMomentum, newPressure);


	/*if(bins[0]->U > 0){
		newDensity[0] = newDensity[0] - deltaT*(densityFlux(0) - densityFlux(rgridNumber-1))/deltaR;
		newMomentum[0] = newMomentum[0] - deltaT*(momentumFlux(0) - momentumFlux(rgridNumber-1) + 0.5*(bins[1]->pressure - bins[rgridNumber-1]->pressure))/deltaR;
		newPressure[0] = newPressure[0]- deltaT*(pressureFlux(0) - pressureFlux(rgridNumber-1) + (gamma-1)*bins[0]->pressure*(volumeFlux(0) - volumeFlux(rgridNumber-1)))/deltaR;
	} else {
		newDensity[0] = newDensity[0] - deltaT*(densityFlux(1) - densityFlux(0))/deltaR;
		newMomentum[0] = newMomentum[0] - deltaT*(momentumFlux(1) - momentumFlux(0) + 0.5*(bins[1]->pressure - bins[0]->pressure))/deltaR;
		newPressure[0] = newPressure[0]- deltaT*(pressureFlux(1) - pressureFlux(0) + (gamma-1)*bins[0]->pressure*(volumeFlux(1) - volumeFlux(0)))/deltaR;
	}
	for(int i = 1; i < rgridNumber - 1; ++i){
		if(bins[i]->U > 0){
			newDensity[i] = newDensity[i] - deltaT*(densityFlux(i) - densityFlux(i-1))/deltaR;
			newMomentum[i] = newMomentum[i] - deltaT*(momentumFlux(i) - momentumFlux(i-1) + (bins[i]->pressure - bins[i-1]->pressure))/deltaR;
			newPressure[i] = newPressure[i]- deltaT*(pressureFlux(i) - pressureFlux(i-1) + (gamma-1)*bins[i]->pressure*(volumeFlux(i) - volumeFlux(i-1)))/deltaR;
		} else {
			newDensity[i] = newDensity[i] - deltaT*(densityFlux(i+1) - densityFlux(i))/deltaR;
			newMomentum[i] = newMomentum[i] - deltaT*(momentumFlux(i+1) - momentumFlux(i) + (bins[i+1]->pressure - bins[i]->pressure))/deltaR;
			newPressure[i] = newPressure[i]- deltaT*(pressureFlux(i+1) - pressureFlux(i) + (gamma-1)*bins[i]->pressure*(volumeFlux(i+1) - volumeFlux(i)))/deltaR;
		}
	}
	if(bins[rgridNumber-1]->U > 0){
		newDensity[rgridNumber-1] = newDensity[rgridNumber-1] - deltaT*(densityFlux(rgridNumber-1) - densityFlux(rgridNumber-2))/deltaR;
		newMomentum[rgridNumber-1] = newMomentum[rgridNumber-1] - deltaT*(momentumFlux(rgridNumber-1) - momentumFlux(rgridNumber-2) + 0.5*(bins[0]->pressure - bins[rgridNumber-2]->pressure))/deltaR;
		newPressure[rgridNumber-1] = newPressure[rgridNumber-1]- deltaT*(pressureFlux(rgridNumber-1) - pressureFlux(rgridNumber-2) + (gamma-1)*bins[rgridNumber-1]->pressure*(volumeFlux(rgridNumber-1) - volumeFlux(rgridNumber-2)))/deltaR;
	} else {
		newDensity[rgridNumber-1] = newDensity[rgridNumber-1] - deltaT*(densityFlux(0) - densityFlux(rgridNumber-1))/deltaR;
		newMomentum[rgridNumber-1] = newMomentum[rgridNumber-1] - deltaT*(momentumFlux(0) - momentumFlux(rgridNumber-1) + 0.5*(bins[0]->pressure - bins[rgridNumber-2]->pressure))/deltaR;
		newPressure[rgridNumber-1] = newPressure[rgridNumber-1]- deltaT*(pressureFlux(0) - pressureFlux(rgridNumber-1) + (gamma-1)*bins[rgridNumber-1]->pressure*(volumeFlux(0) - volumeFlux(rgridNumber-1)))/deltaR;
	}*/


	double* densityFluxes = new double[rgridNumber];
	double* momentumFluxes = new double[rgridNumber];
	double* pressureFluxes = new double[rgridNumber];

	double* momentumPressureFluxes = new double[rgridNumber];
	double* volumeFluxes = new double[rgridNumber];

	for(int i = 0; i < rgridNumber; ++i){
		newDensity[i] = bins[i]->density;
		newMomentum[i] = bins[i]->density*bins[i]->U;
		newPressure[i] = bins[i]->pressure;

		/*tempDensity[i] = bins[i]->density;
		tempMomentum[i] = bins[i]->density*bins[i]->U;
		tempPressure[i] = bins[i]->pressure;*/

		densityFluxes[i] = densityFlux(i);
		momentumFluxes[i] = momentumFlux(i);
		pressureFluxes[i] = pressureFlux(i);
		momentumPressureFluxes[i] = momentumPressureFlux(i);
		volumeFluxes[i] = (gamma-1)*volumeFlux(i);
	}
	//runge-kutt step 
	/*deltaT = deltaT/2;

	tvd(tempDensity, densityFluxes, c);
	tvd(tempMomentum, momentumFluxes, c);
	tvd(tempPressure, pressureFluxes, c);

	for(int i = 0; i < rgridNumber; ++i){
		tempDensity[i] = tempDensity[i];
		tempMomentum[i] = tempMomentum[i];
		if((tempPressure[i] != tempPressure[i]) || (0*tempPressure[i] != 0*tempPressure[i])){
			printf("NaN pressure\n");
			Sleep(500);
		}
		tempPressure[i] = tempPressure[i]/(bins[i]->pressure);
	}
	tvd(tempMomentum, momentumPressureFluxes, c);
	tvd(tempPressure, volumeFluxes, c);
	for(int i = 0; i < rgridNumber; ++i){
		tempMomentum[i] = tempMomentum[i];
		if((tempPressure[i] != tempPressure[i]) || (0*tempPressure[i] != 0*tempPressure[i])){
			printf("NaN pressure\n");
			Sleep(500);
		}
		if((bins[i]->pressure*tempPressure[i] != tempPressure[i]*bins[i]->pressure) || (0*tempPressure[i]*bins[i]->pressure != 0*tempPressure[i]*bins[i]->pressure)){
			printf("NaN pressure\n");
			Sleep(500);
		}
		tempPressure[i] = tempPressure[i]*bins[i]->pressure;
	}

	for(int i = 0; i < rgridNumber; ++i){
		if((tempMomentum[i] != tempMomentum[i]) || (0*tempMomentum[i] != 0*tempMomentum[i])){
			printf("NaN velocity\n");
			Sleep(500);
		}
		bins[i]->U = (tempMomentum[i]/bins[i]->density);
		if(tempDensity[i] < 0){
			bins[i]->density = epsilon;
			printf("density < 0\n");
		} else {
			if((tempDensity[i] != tempDensity[i]) || (0*tempDensity[i] != 0*tempDensity[i])){
				printf("NaN density\n");
				Sleep(1000);
			}
			bins[i]->density = tempDensity[i];
		}
		if(bins[i]->density < 100*epsilon*density0){
			bins[i]->U = 0;
		}
		bins[i]->pressure = tempPressure[i];
		if(bins[i]->pressure < 0){
			bins[i]->pressure = epsilon;
			printf("pressure < 0\n");
		} else {
			if((bins[i]->pressure != bins[i]->pressure) || (0*bins[i]->pressure != 0*bins[i]->pressure)){
				printf("NaN pressure\n");
				Sleep(500);
			}
		}
		if(bins[i]->density < 100*epsilon*density0){
			bins[i]->pressure = epsilon*tempPressure[i];
		}
	}
	deltaT = deltaT*2;
	// end of runge kutt
	*/

	for(int i = 0; i < rgridNumber; ++i){
		densityFluxes[i] = densityFlux(i);
		momentumFluxes[i] = momentumFlux(i);
		pressureFluxes[i] = pressureFlux(i);
		momentumPressureFluxes[i] = momentumPressureFlux(i);
		volumeFluxes[i] = (gamma-1)*volumeFlux(i);
	}

	convectionTVD(newDensity, densityFluxes);
	convectionTVD(newMomentum, momentumFluxes);
	convectionTVD(newPressure, pressureFluxes);

	for(int i = 0; i < rgridNumber; ++i){
		newDensity[i] = newDensity[i];
		newMomentum[i] = newMomentum[i];
		newPressure[i] = newPressure[i]/(bins[i]->pressure);
	}
	tvd(newMomentum, momentumPressureFluxes, c);
	tvd(newPressure, volumeFluxes, c);
	for(int i = 0; i < rgridNumber; ++i){
		newMomentum[i] = newMomentum[i];
		newPressure[i] = newPressure[i]*bins[i]->pressure;
	}

	for(int i = 0; i < rgridNumber - 1; ++i){
		tempDensity[i] = (1 - scaleParameter)*newDensity[i] + scaleParameter*newDensity[i+1];
		tempMomentum[i] = (1 - scaleParameter)*newMomentum[i] + scaleParameter*newMomentum[i+1];
		tempPressure[i] = (1 - scaleParameter)*newPressure[i] + scaleParameter*newPressure[i+1];
	}
	tempDensity[rgridNumber - 1] = (1 - scaleParameter)*newDensity[rgridNumber-1] + scaleParameter*newDensity[0];
	tempMomentum[rgridNumber - 1] = (1 - scaleParameter)*newMomentum[rgridNumber-1] + scaleParameter*newMomentum[0];
	tempPressure[rgridNumber - 1] = (1 - scaleParameter)*newPressure[rgridNumber-1] + scaleParameter*newPressure[0];

	for(int i = 0; i < rgridNumber; ++i){
		newDensity[i] = tempDensity[i];
		newMomentum[i] = tempMomentum[i] - newDensity[i]*c;
		newPressure[i] = tempPressure[i];
	}




	delete[] densityFluxes;
	delete[] momentumFluxes;
	delete[] pressureFluxes;

	delete[] tempDensity;
	delete[] tempMomentum;
	delete[] tempPressure;

	delete[] momentumPressureFluxes;
	delete[] volumeFluxes;
}

double Simulation::densityFlux(int i){
	if(i == -1){
		return 0;
	}
	double flux = bins[i]->U*bins[i]->density;
	if(flux != flux || (0*flux != 0*flux)){
		printf("NaN densityFlux\n");
	}
	return flux;
}

double Simulation::momentumFlux(int i){
	if(i == -1){
		return 0;
	}
	double flux = bins[i]->density*bins[i]->U*bins[i]->U;
	if(flux != flux || (0*flux != 0*flux)){
		printf("NaN momentumFlux\n");
	}
	return flux;
}

double Simulation::pressureFlux(int i){
	if(i == -1){
		return 0;
	}
	double flux = bins[i]->U*bins[i]->pressure;
	if(flux != flux || (0*flux != 0*flux)){
		printf("NaN energyFlux\n");
	}
	return flux;
}

double Simulation::momentumPressureFlux(int i){
	if(i == -1){
		return 0;
	}
	if(bins[i]->pressure != bins[i]->pressure || (0*bins[i]->pressure != 0*bins[i]->pressure)){
		printf("NaN pressureFlux\n");
	}
	return bins[i]->pressure;
}

double Simulation::volumeFlux(int i){
	double flux = bins[i]->U;
	if(flux != flux || (0*flux != 0*flux)){
		printf("NaN volumeFlux\n");
	}
	return flux;
}

double Simulation::fullMomentumFlux(int i){
	return (bins[i]->density*bins[i]->U*bins[i]->U + bins[i]->pressure);
}

double Simulation::fullEnergyFlux(int i){
	return (bins[i]->U*(bins[i]->pressure/(gamma-1) + bins[i]->density*bins[i]->U*bins[i]->U/2 + bins[i]->pressure));
}

double Simulation::vanleer(double a, double b){
	if(abs(a + b) < epsilon){
		return 0;
	} else if( a*b < 0) {
		return 0;
	} else {
		return 2*a*b/(a + b);
	}
}

double Simulation::findMaxVelocity(){
	double maxVelocity = 0;
	for(int i = 0; i < rgridNumber; ++i){
		if( abs(bins[i]->U) + bins[i]->soundSpeed() > abs(maxVelocity)){
			maxVelocity = abs(bins[i]->U) + bins[i]->soundSpeed();
		}
	}
	return maxVelocity;
}

void Simulation::tvd(double* value, double* valueFlux, double maxVelocity){
	double* fr = new double[rgridNumber];
	double* fl = new double[rgridNumber];
    double* flux = new double[rgridNumber];

	for(int i = 0; i < rgridNumber - 1; ++i){
		fr[i] = (valueFlux[i] + valueFlux[i+1])/2;
		fl[i+1] = fr[i];
	}
	fr[rgridNumber - 1] = (valueFlux[rgridNumber - 1] + valueFlux[0])/2;
	fl[0] = fr[rgridNumber - 1];

	/*for(int i = 0; i < rgridNumber - 1; ++i){
		if(bins[i]->U > 0){
			if(bins[i+1]->U > 0){
				fr[i] = valueFlux[i];
				fl[i+1] = fr[i];
			} else {
				fr[i] = (valueFlux[i] + valueFlux[i+1])/2;
				fl[i+1] = fr[i];
			}
		} else {
			if(bins[i+1]->U > 0){
				fr[i] = (valueFlux[i] + valueFlux[i+1])/2;
				fl[i+1] = fr[i];
			} else {				
				fr[i] = valueFlux[i+1];
				fl[i+1] = fr[i];
			}
		}

		
	}

	if(bins[rgridNumber - 1]->U > 0){
		if(bins[0]->U > 0){
			fr[rgridNumber - 1] = valueFlux[rgridNumber - 1];
			fl[0] = fr[rgridNumber - 1];
		} else {
			fr[rgridNumber - 1] = (valueFlux[rgridNumber - 1] + valueFlux[0])/2;
			fl[0] = fr[rgridNumber - 1];
		}
	} else {
		if(bins[0]->U > 0){
			fr[rgridNumber - 1] = (valueFlux[rgridNumber - 1] + valueFlux[0])/2;
			fl[0] = fr[rgridNumber - 1];
		} else {				
			fr[rgridNumber - 1] = valueFlux[0];
			fl[0] = fr[rgridNumber - 1];
		}
	}*/

	for(int i = 0; i < rgridNumber; ++i){
		value[i] = value[i] - deltaT*(fr[i] - fl[i])/deltaR;
		if(value[i] != value[i] || (0*value[i] != 0*value[i])){
		    printf("NaN value");
		}
	}

	/*double* tempValue = new double[rgridNumber];

	if(maxVelocity > 0){
	  fr[0] = maxVelocity*value[0] + valueFlux[0];
	  for(int i = 1; i < rgridNumber; ++i){
		  fr[i] = maxVelocity*value[i] + valueFlux[i];
		  fl[i] = maxVelocity*value[i-1] - valueFlux[i-1];
	  }
      fl[0] = maxVelocity*value[rgridNumber - 1] - valueFlux[rgridNumber - 1];

	  for(int i = 0; i < rgridNumber; ++i){
		  flux[i] = (fr[i] - fl[i])/2;
	  }

	  tempValue[0] = value[0] - deltaT*0.5*(flux[0] - flux[rgridNumber - 1])/deltaR;
	  for(int i = 1; i < rgridNumber; ++i){
		  tempValue[i] = value[i] - deltaT*0.5*(flux[i] - flux[i-1])/deltaR;
	  }

	  fr[0] = maxVelocity*tempValue[0] + valueFlux[0];
	  for(int i = 1; i < rgridNumber; ++i){
		  fr[i] = maxVelocity*tempValue[i] + valueFlux[i];
		  fl[i] = maxVelocity*tempValue[i-1] - valueFlux[i-1];
	  }
      fl[0] = maxVelocity*tempValue[rgridNumber - 1] - valueFlux[rgridNumber - 1];

	  flux[0] = (fr[0] - fl[0])/2;
	  for(int i = 1; i < rgridNumber - 1; ++i){
		  /*double dfrp = (fr[i+1] - fr[i])/2;
		  double dfrm = (fr[i] - fr[i-1])/2;
		  double dfr = vanleer(dfrp, dfrm);

		  double dflp = (fl[i]-fl[i+1])/2;
		  double dflm = (fl[i-1] - fl[i])/2;
		  double dfl = vanleer(dflp, dflm);

		  flux[i] = (fr[i] - fl[i] + dfr - dfl)/2;*/
		  /*flux[i] = (fr[i] - fl[i])/2;
	  }
	  flux[rgridNumber - 1] = (fr[rgridNumber - 1] - fl[rgridNumber - 1])/2;


	  value[0] = value[0] - deltaT*(flux[0] - flux[rgridNumber - 1])/deltaR;
	  if(value[0] != value[0] || (0*value[0] != 0*value[0])){
	    printf("NaN value");
	  }
	  for(int i = 1; i < rgridNumber; ++i){
		  value[i] = value[i] - deltaT*(flux[i] - flux[i-1])/deltaR;
		  if(value[i] != value[i] || (0*value[i] != 0*value[i])){
			  printf("NaN value");
		  }
	  }
	} else {
	  for(int i = 0; i < rgridNumber - 1; ++i){
		  fr[i] = maxVelocity*value[i+1] + valueFlux[i+1];
		  fl[i] = maxVelocity*value[i] - valueFlux[i];
	  }
	  fr[rgridNumber - 1] = maxVelocity*value[0] + valueFlux[0];
	  fl[rgridNumber - 1] = maxVelocity*value[rgridNumber - 1] - valueFlux[rgridNumber - 1];

	  for(int i = 0; i < rgridNumber; ++i){
		  flux[i] = (fr[i] - fl[i])/2;
	  }

	  
	  for(int i = 0; i < rgridNumber - 1; ++i){
		  tempValue[i] = value[i] - deltaT*0.5*(flux[i+1] - flux[i])/deltaR;
	  }
	  tempValue[rgridNumber - 1] = value[rgridNumber - 1] - deltaT*0.5*(flux[0] - flux[rgridNumber - 1])/deltaR;

	  for(int i = 0; i < rgridNumber-1; ++i){
		  fr[i] = maxVelocity*tempValue[i+1] + valueFlux[i+1];
		  fl[i] = maxVelocity*tempValue[i] - valueFlux[i];
	  }
	  fr[rgridNumber- 1] = maxVelocity*tempValue[0] + valueFlux[0];
	  fl[rgridNumber- 1] = maxVelocity*tempValue[rgridNumber - 1] - valueFlux[rgridNumber - 1];

	  flux[0] = (fr[0] - fl[0])/2;
	  for(int i = 0; i < rgridNumber - 2; ++i){
	 	  /*double dfrp = -(fr[i+2] - fr[i+1])/2;
		  double dfrm = -(fr[i+1] - fr[i])/2;
		  double dfr = vanleer(dfrp, dfrm);

		  double dflp = (fl[i+2] - fl[i+1])/2;
		  double dflm = (fl[i+1] - fl[i])/2;
		  double dfl = vanleer(dflp, dflm);
		  
		  flux[i] = (fr[i] - fl[i] + dfr - dfl)/2;*/
		  /*flux[i] = (fr[i] - fl[i])/2;
	  }
	  flux[rgridNumber - 2] = (fr[rgridNumber - 2] - fl[rgridNumber - 2])/2;
	  flux[rgridNumber - 1] = (fr[rgridNumber - 1] - fl[rgridNumber - 1])/2;

	  for(int i = 0; i < rgridNumber - 1; ++i){
		  value[i] = value[i] - deltaT*(flux[i+1] - flux[i])/deltaR;
		  if(value[i] != value[i] || (0*value[i] != 0*value[i])){
			  printf("NaN value");
		  }
	  }
	  value[rgridNumber - 1] = value[rgridNumber - 1] - deltaT*(flux[0] - flux[rgridNumber - 1])/deltaR;
	}*/

	delete[] fr;
	delete[] fl;
	delete[] flux;
	//delete[] tempValue;
}

void Simulation::convectionTVD(double* value, double* valueFlux){
	double* fr = new double[rgridNumber];
	double* fl = new double[rgridNumber];
    double* flux = new double[rgridNumber];

	for(int i = 0; i < rgridNumber - 1; ++i){
		if(bins[i]->U > 0){
			if(bins[i+1]->U > 0){
				fr[i] = valueFlux[i];
				fl[i+1] = fr[i];
			} else {
				fr[i] = (valueFlux[i] + valueFlux[i+1])/2;
				fl[i+1] = fr[i];
			}
		} else {
			if(bins[i+1]->U > 0){
				fr[i] = (valueFlux[i] + valueFlux[i+1])/2;;
				fl[i+1] = fr[i];
			} else {				
				fr[i] = valueFlux[i+1];
				fl[i+1] = fr[i];
			}
		}
	}

	if(bins[rgridNumber - 1]->U > 0){
		if(bins[0]->U > 0){
			fr[rgridNumber - 1] = valueFlux[rgridNumber - 1];
			fl[0] = fr[rgridNumber - 1];
		} else {
			fr[rgridNumber - 1] = (valueFlux[rgridNumber - 1] + valueFlux[0])/2;
			fl[0] = fr[rgridNumber - 1];
		}
	} else {
		if(bins[0]->U > 0){
			fr[rgridNumber - 1] = (valueFlux[rgridNumber - 1] + valueFlux[0])/2;
			fl[0] = fr[rgridNumber - 1];
		} else {				
			fr[rgridNumber - 1] = valueFlux[0];
			fl[0] = fr[rgridNumber - 1];
		}
	}

	for(int i = 0; i < rgridNumber; ++i){
		value[i] = value[i] - deltaT*(fr[i] - fl[i])/deltaR;
		if(value[i] != value[i] || (0*value[i] != 0*value[i])){
		    printf("NaN value");
		}
	}

	delete[] fr;
	delete[] fl;
	delete[] flux;
}

double Simulation::kurzrockMinDeltaT(){
	double minDeltaT = abs(bins[0]->U)*deltaR/((abs(bins[0]->U) + sqrt(2.0)*bins[0]->soundSpeed())*(abs(bins[0]->U) + sqrt(2.0)*bins[0]->soundSpeed()));
	for(int i = 1; i < rgridNumber; ++i){
		if(abs(bins[i]->U)*deltaR/((abs(bins[i]->U) + sqrt(2.0)*bins[i]->soundSpeed())*(abs(bins[i]->U) + sqrt(2.0)*bins[i]->soundSpeed())) < minDeltaT){
			minDeltaT = abs(bins[i]->U)*deltaR/((abs(bins[i]->U) + sqrt(2.0)*bins[i]->soundSpeed())*(abs(bins[i]->U) + sqrt(2.0)*bins[i]->soundSpeed()));
		}
	}
	return minDeltaT;
}

void Simulation::LaxVendorf(double* newDensity, double* newMomentum, double* newEnergy){
	// P.Rouch Computative Hydrodyamic 5.5.5

	//modifization for u > < 0

		/*newDensity[0] = newDensity[0] - deltaT*0.5*(densityFlux(1) - densityFlux(rgridNumber-1))/deltaR 
			+ 0.25*deltaT*deltaT*((2*(fullMomentumFlux(1) - fullMomentumFlux(0)) - 2*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1))))/(deltaR*deltaR);
		newMomentum[0] = newMomentum[0] - deltaT*0.5*(fullMomentumFlux(1) - fullMomentumFlux(rgridNumber-1))/deltaR 
			+ 0.25*deltaT*deltaT*((-0.5*(3-gamma)*(bins[1]->U*bins[1]->U + bins[0]->U*bins[0]->U)*(densityFlux(1) - densityFlux(0)) + 0.5*(3-gamma)*(bins[0]->U*bins[0]->U + bins[rgridNumber-1]->U*bins[rgridNumber-1]->U)*(densityFlux(0) - densityFlux(rgridNumber-1))) 
			+ (2*(3-gamma)*(bins[1]->U + bins[0]->U)*(fullMomentumFlux(1) - fullMomentumFlux(0)) - 2*(3-gamma)*(bins[0]->U+bins[rgridNumber-1]->U)*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1))) 
			+ (2*(gamma - 1)*(fullEnergyFlux(1) - fullEnergyFlux(0)) - 2*(gamma - 1)*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-1))))/(deltaR*deltaR);
		newEnergy[0] = newEnergy[0] - deltaT*0.5*(fullEnergyFlux(1) - fullEnergyFlux(rgridNumber-1))/deltaR 
			+ 0.25*deltaT*deltaT*(((-gamma*bins[1]->U*bins[1]->fullEnergy()/bins[1]->density + (gamma - 1)*bins[0+1]->U*bins[0+1]->U*bins[0+1]->U - gamma*bins[0]->U*bins[0]->fullEnergy()/bins[0]->density + (gamma - 1)*bins[0]->U*bins[0]->U*bins[0]->U)*(densityFlux(0+1) - densityFlux(0)) - (-gamma*bins[0]->U*bins[0]->fullEnergy()/bins[0]->density + (gamma - 1)*bins[0]->U*bins[0]->U*bins[0]->U - gamma*bins[rgridNumber-1]->U*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density + (gamma - 1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U)*(densityFlux(0) - densityFlux(rgridNumber-1))) 
			+ ((gamma*bins[1]->fullEnergy()/bins[1]->density - 3*(gamma-1)*bins[1]->U*bins[1]->U/2 + gamma*bins[0]->fullEnergy()/bins[0]->density - 3*(gamma-1)*bins[0]->U*bins[0]->U/2)*(fullMomentumFlux(1) - fullMomentumFlux(0)) - (gamma*bins[0]->fullEnergy()/bins[0]->density - 3*(gamma-1)*bins[0]->U*bins[0]->U/2 + gamma*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density - 3*(gamma-1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U/2)*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1))) 
			+ (gamma*(bins[1]->U + bins[0]->U)*(fullEnergyFlux(1) - fullEnergyFlux(0)) - gamma*(bins[0]->U + bins[rgridNumber-1]->U)*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-1))))/(deltaR*deltaR);*/
	if(bins[0]->U >= 0){
		newDensity[0] = newDensity[0] - deltaT*(densityFlux(0) - densityFlux(rgridNumber-1))/deltaR 
			+ 0.25*deltaT*deltaT*((2*(fullMomentumFlux(1) - fullMomentumFlux(0)) - 2*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1))))/(deltaR*deltaR);
		newMomentum[0] = newMomentum[0] - deltaT*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1))/deltaR 
			+ 0.25*deltaT*deltaT*((-0.5*(3-gamma)*(bins[1]->U*bins[1]->U + bins[0]->U*bins[0]->U)*(densityFlux(1) - densityFlux(0)) + 0.5*(3-gamma)*(bins[0]->U*bins[0]->U + bins[rgridNumber-1]->U*bins[rgridNumber-1]->U)*(densityFlux(0) - densityFlux(rgridNumber-1))) 
			+ (2*(3-gamma)*(bins[1]->U + bins[0]->U)*(fullMomentumFlux(1) - fullMomentumFlux(0)) - 2*(3-gamma)*(bins[0]->U+bins[rgridNumber-1]->U)*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1))) 
			+ (2*(gamma - 1)*(fullEnergyFlux(1) - fullEnergyFlux(0)) - 2*(gamma - 1)*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-1))))/(deltaR*deltaR);
		newEnergy[0] = newEnergy[0] - deltaT*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-1))/deltaR 
			+ 0.25*deltaT*deltaT*(((-gamma*bins[1]->U*bins[1]->fullEnergy()/bins[1]->density + (gamma - 1)*bins[0+1]->U*bins[0+1]->U*bins[0+1]->U - gamma*bins[0]->U*bins[0]->fullEnergy()/bins[0]->density + (gamma - 1)*bins[0]->U*bins[0]->U*bins[0]->U)*(densityFlux(0+1) - densityFlux(0)) - (-gamma*bins[0]->U*bins[0]->fullEnergy()/bins[0]->density + (gamma - 1)*bins[0]->U*bins[0]->U*bins[0]->U - gamma*bins[rgridNumber-1]->U*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density + (gamma - 1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U)*(densityFlux(0) - densityFlux(rgridNumber-1))) 
			+ ((gamma*bins[1]->fullEnergy()/bins[1]->density - 3*(gamma-1)*bins[1]->U*bins[1]->U/2 + gamma*bins[0]->fullEnergy()/bins[0]->density - 3*(gamma-1)*bins[0]->U*bins[0]->U/2)*(fullMomentumFlux(1) - fullMomentumFlux(0)) - (gamma*bins[0]->fullEnergy()/bins[0]->density - 3*(gamma-1)*bins[0]->U*bins[0]->U/2 + gamma*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density - 3*(gamma-1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U/2)*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1))) 
			+ (gamma*(bins[1]->U + bins[0]->U)*(fullEnergyFlux(1) - fullEnergyFlux(0)) - gamma*(bins[0]->U + bins[rgridNumber-1]->U)*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-1))))/(deltaR*deltaR);
	} else {
		newDensity[0] = newDensity[0] - deltaT*(densityFlux(1) - densityFlux(0))/deltaR 
			+ 0.25*deltaT*deltaT*((2*(fullMomentumFlux(1) - fullMomentumFlux(0)) - 2*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1))))/(deltaR*deltaR);
		newMomentum[0] = newMomentum[0] - deltaT*(fullMomentumFlux(1) - fullMomentumFlux(0))/deltaR 
			+ 0.25*deltaT*deltaT*((-0.5*(3-gamma)*(bins[1]->U*bins[1]->U + bins[0]->U*bins[0]->U)*(densityFlux(1) - densityFlux(0)) + 0.5*(3-gamma)*(bins[0]->U*bins[0]->U + bins[rgridNumber-1]->U*bins[rgridNumber-1]->U)*(densityFlux(0) - densityFlux(rgridNumber-1))) 
			+ (2*(3-gamma)*(bins[1]->U + bins[0]->U)*(fullMomentumFlux(1) - fullMomentumFlux(0)) - 2*(3-gamma)*(bins[0]->U+bins[rgridNumber-1]->U)*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1))) 
			+ (2*(gamma - 1)*(fullEnergyFlux(1) - fullEnergyFlux(0)) - 2*(gamma - 1)*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-1))))/(deltaR*deltaR);
		newEnergy[0] = newEnergy[0] - deltaT*(fullEnergyFlux(1) - fullEnergyFlux(0))/deltaR 
			+ 0.25*deltaT*deltaT*(((-gamma*bins[1]->U*bins[1]->fullEnergy()/bins[1]->density + (gamma - 1)*bins[0+1]->U*bins[0+1]->U*bins[0+1]->U - gamma*bins[0]->U*bins[0]->fullEnergy()/bins[0]->density + (gamma - 1)*bins[0]->U*bins[0]->U*bins[0]->U)*(densityFlux(0+1) - densityFlux(0)) - (-gamma*bins[0]->U*bins[0]->fullEnergy()/bins[0]->density + (gamma - 1)*bins[0]->U*bins[0]->U*bins[0]->U - gamma*bins[rgridNumber-1]->U*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density + (gamma - 1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U)*(densityFlux(0) - densityFlux(rgridNumber-1))) 
			+ ((gamma*bins[1]->fullEnergy()/bins[1]->density - 3*(gamma-1)*bins[1]->U*bins[1]->U/2 + gamma*bins[0]->fullEnergy()/bins[0]->density - 3*(gamma-1)*bins[0]->U*bins[0]->U/2)*(fullMomentumFlux(1) - fullMomentumFlux(0)) - (gamma*bins[0]->fullEnergy()/bins[0]->density - 3*(gamma-1)*bins[0]->U*bins[0]->U/2 + gamma*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density - 3*(gamma-1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U/2)*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1))) 
			+ (gamma*(bins[1]->U + bins[0]->U)*(fullEnergyFlux(1) - fullEnergyFlux(0)) - gamma*(bins[0]->U + bins[rgridNumber-1]->U)*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-1))))/(deltaR*deltaR);
	}
			
	for(int i = 1; i < rgridNumber - 1; ++i){


				/*newDensity[i] = newDensity[i] - deltaT*0.5*(densityFlux(i+1) - densityFlux(i-1))/deltaR 
			+ 0.25*deltaT*deltaT*((2*(fullMomentumFlux(i+1) - fullMomentumFlux(i)) - 2*(fullMomentumFlux(i) - fullMomentumFlux(i-1))))/(deltaR*deltaR);
		newMomentum[i] = newMomentum[i] - deltaT*0.5*(fullMomentumFlux(i+1) - fullMomentumFlux(i-1))/deltaR 
			+ 0.25*deltaT*deltaT*((-0.5*(3-gamma)*(bins[i+1]->U*bins[i+1]->U + bins[i]->U*bins[i]->U)*(densityFlux(i+1) - densityFlux(i)) + 0.5*(3-gamma)*(bins[i]->U*bins[i]->U + bins[i-1]->U*bins[i-1]->U)*(densityFlux(i) - densityFlux(i-1))) 
			+ (2*(3-gamma)*(bins[i+1]->U + bins[i]->U)*(fullMomentumFlux(i+1) - fullMomentumFlux(i)) - 2*(3-gamma)*(bins[i]->U+bins[i-1]->U)*(fullMomentumFlux(i) - fullMomentumFlux(i-1))) 
			+ (2*(gamma - 1)*(fullEnergyFlux(i+1) - fullEnergyFlux(i)) - 2*(gamma - 1)*(fullEnergyFlux(i) - fullEnergyFlux(i-1))))/(deltaR*deltaR);
		newEnergy[i] = newEnergy[i] - deltaT*0.5*(fullEnergyFlux(i+1) - fullEnergyFlux(i-1))/deltaR 
			+ 0.25*deltaT*deltaT*(((-gamma*bins[i+1]->U*bins[i+1]->fullEnergy()/bins[i+1]->density + (gamma - 1)*bins[i+1]->U*bins[i+1]->U*bins[i+1]->U - gamma*bins[i]->U*bins[i]->fullEnergy()/bins[i]->density + (gamma - 1)*bins[i]->U*bins[i]->U*bins[i]->U)*(densityFlux(i+1) - densityFlux(i)) - (-gamma*bins[i]->U*bins[i]->fullEnergy()/bins[i]->density + (gamma - 1)*bins[i]->U*bins[i]->U*bins[i]->U - gamma*bins[i-1]->U*bins[i-1]->fullEnergy()/bins[i-1]->density + (gamma - 1)*bins[i-1]->U*bins[i-1]->U*bins[i-1]->U)*(densityFlux(i) - densityFlux(i-1))) 
			+ ((gamma*bins[i+1]->fullEnergy()/bins[i+1]->density - 3*(gamma-1)*bins[i+1]->U*bins[i+1]->U/2 + gamma*bins[i]->fullEnergy()/bins[i]->density - 3*(gamma-1)*bins[i]->U*bins[i]->U/2)*(fullMomentumFlux(i+1) - fullMomentumFlux(i)) - (gamma*bins[i]->fullEnergy()/bins[i]->density - 3*(gamma-1)*bins[i]->U*bins[i]->U/2 + gamma*bins[i-1]->fullEnergy()/bins[i-1]->density - 3*(gamma-1)*bins[i-1]->U*bins[i-1]->U/2)*(fullMomentumFlux(i) - fullMomentumFlux(i-1))) 
			+ (gamma*(bins[i+1]->U + bins[i]->U)*(fullEnergyFlux(i+1) - fullEnergyFlux(i)) - gamma*(bins[i]->U + bins[i-1]->U)*(fullEnergyFlux(i) - fullEnergyFlux(i-1))))/(deltaR*deltaR);*/
		if(bins[i]->U >= 0){
		newDensity[i] = newDensity[i] - deltaT*(densityFlux(i) - densityFlux(i-1))/deltaR 
			+ 0.25*deltaT*deltaT*((2*(fullMomentumFlux(i+1) - fullMomentumFlux(i)) - 2*(fullMomentumFlux(i) - fullMomentumFlux(i-1))))/(deltaR*deltaR);
		newMomentum[i] = newMomentum[i] - deltaT*(fullMomentumFlux(i) - fullMomentumFlux(i-1))/deltaR 
			+ 0.25*deltaT*deltaT*((-0.5*(3-gamma)*(bins[i+1]->U*bins[i+1]->U + bins[i]->U*bins[i]->U)*(densityFlux(i+1) - densityFlux(i)) + 0.5*(3-gamma)*(bins[i]->U*bins[i]->U + bins[i-1]->U*bins[i-1]->U)*(densityFlux(i) - densityFlux(i-1))) 
			+ (2*(3-gamma)*(bins[i+1]->U + bins[i]->U)*(fullMomentumFlux(i+1) - fullMomentumFlux(i)) - 2*(3-gamma)*(bins[i]->U+bins[i-1]->U)*(fullMomentumFlux(i) - fullMomentumFlux(i-1))) 
			+ (2*(gamma - 1)*(fullEnergyFlux(i+1) - fullEnergyFlux(i)) - 2*(gamma - 1)*(fullEnergyFlux(i) - fullEnergyFlux(i-1))))/(deltaR*deltaR);
		newEnergy[i] = newEnergy[i] - deltaT*(fullEnergyFlux(i) - fullEnergyFlux(i-1))/deltaR 
			+ 0.25*deltaT*deltaT*(((-gamma*bins[i+1]->U*bins[i+1]->fullEnergy()/bins[i+1]->density + (gamma - 1)*bins[i+1]->U*bins[i+1]->U*bins[i+1]->U - gamma*bins[i]->U*bins[i]->fullEnergy()/bins[i]->density + (gamma - 1)*bins[i]->U*bins[i]->U*bins[i]->U)*(densityFlux(i+1) - densityFlux(i)) - (-gamma*bins[i]->U*bins[i]->fullEnergy()/bins[i]->density + (gamma - 1)*bins[i]->U*bins[i]->U*bins[i]->U - gamma*bins[i-1]->U*bins[i-1]->fullEnergy()/bins[i-1]->density + (gamma - 1)*bins[i-1]->U*bins[i-1]->U*bins[i-1]->U)*(densityFlux(i) - densityFlux(i-1))) 
			+ ((gamma*bins[i+1]->fullEnergy()/bins[i+1]->density - 3*(gamma-1)*bins[i+1]->U*bins[i+1]->U/2 + gamma*bins[i]->fullEnergy()/bins[i]->density - 3*(gamma-1)*bins[i]->U*bins[i]->U/2)*(fullMomentumFlux(i+1) - fullMomentumFlux(i)) - (gamma*bins[i]->fullEnergy()/bins[i]->density - 3*(gamma-1)*bins[i]->U*bins[i]->U/2 + gamma*bins[i-1]->fullEnergy()/bins[i-1]->density - 3*(gamma-1)*bins[i-1]->U*bins[i-1]->U/2)*(fullMomentumFlux(i) - fullMomentumFlux(i-1))) 
			+ (gamma*(bins[i+1]->U + bins[i]->U)*(fullEnergyFlux(i+1) - fullEnergyFlux(i)) - gamma*(bins[i]->U + bins[i-1]->U)*(fullEnergyFlux(i) - fullEnergyFlux(i-1))))/(deltaR*deltaR);
		} else {
		newDensity[i] = newDensity[i] - deltaT*(densityFlux(i+1) - densityFlux(i))/deltaR 
			+ 0.25*deltaT*deltaT*((2*(fullMomentumFlux(i+1) - fullMomentumFlux(i)) - 2*(fullMomentumFlux(i) - fullMomentumFlux(i-1))))/(deltaR*deltaR);
		newMomentum[i] = newMomentum[i] - deltaT*(fullMomentumFlux(i+1) - fullMomentumFlux(i))/deltaR 
			+ 0.25*deltaT*deltaT*((-0.5*(3-gamma)*(bins[i+1]->U*bins[i+1]->U + bins[i]->U*bins[i]->U)*(densityFlux(i+1) - densityFlux(i)) + 0.5*(3-gamma)*(bins[i]->U*bins[i]->U + bins[i-1]->U*bins[i-1]->U)*(densityFlux(i) - densityFlux(i-1))) 
			+ (2*(3-gamma)*(bins[i+1]->U + bins[i]->U)*(fullMomentumFlux(i+1) - fullMomentumFlux(i)) - 2*(3-gamma)*(bins[i]->U+bins[i-1]->U)*(fullMomentumFlux(i) - fullMomentumFlux(i-1))) 
			+ (2*(gamma - 1)*(fullEnergyFlux(i+1) - fullEnergyFlux(i)) - 2*(gamma - 1)*(fullEnergyFlux(i) - fullEnergyFlux(i-1))))/(deltaR*deltaR);
		newEnergy[i] = newEnergy[i] - deltaT*(fullEnergyFlux(i+1) - fullEnergyFlux(i))/deltaR 
			+ 0.25*deltaT*deltaT*(((-gamma*bins[i+1]->U*bins[i+1]->fullEnergy()/bins[i+1]->density + (gamma - 1)*bins[i+1]->U*bins[i+1]->U*bins[i+1]->U - gamma*bins[i]->U*bins[i]->fullEnergy()/bins[i]->density + (gamma - 1)*bins[i]->U*bins[i]->U*bins[i]->U)*(densityFlux(i+1) - densityFlux(i)) - (-gamma*bins[i]->U*bins[i]->fullEnergy()/bins[i]->density + (gamma - 1)*bins[i]->U*bins[i]->U*bins[i]->U - gamma*bins[i-1]->U*bins[i-1]->fullEnergy()/bins[i-1]->density + (gamma - 1)*bins[i-1]->U*bins[i-1]->U*bins[i-1]->U)*(densityFlux(i) - densityFlux(i-1))) 
			+ ((gamma*bins[i+1]->fullEnergy()/bins[i+1]->density - 3*(gamma-1)*bins[i+1]->U*bins[i+1]->U/2 + gamma*bins[i]->fullEnergy()/bins[i]->density - 3*(gamma-1)*bins[i]->U*bins[i]->U/2)*(fullMomentumFlux(i+1) - fullMomentumFlux(i)) - (gamma*bins[i]->fullEnergy()/bins[i]->density - 3*(gamma-1)*bins[i]->U*bins[i]->U/2 + gamma*bins[i-1]->fullEnergy()/bins[i-1]->density - 3*(gamma-1)*bins[i-1]->U*bins[i-1]->U/2)*(fullMomentumFlux(i) - fullMomentumFlux(i-1))) 
			+ (gamma*(bins[i+1]->U + bins[i]->U)*(fullEnergyFlux(i+1) - fullEnergyFlux(i)) - gamma*(bins[i]->U + bins[i-1]->U)*(fullEnergyFlux(i) - fullEnergyFlux(i-1))))/(deltaR*deltaR);
		}

	}

		/*newDensity[rgridNumber-1] = newDensity[rgridNumber-1] - deltaT*0.5*(densityFlux(0) - densityFlux(rgridNumber-2))/deltaR 
			+ 0.25*deltaT*deltaT*((2*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1)) - 2*(fullMomentumFlux(rgridNumber-1) - fullMomentumFlux(rgridNumber-2))))/(deltaR*deltaR);
		newMomentum[rgridNumber-1] = newMomentum[rgridNumber-1] - deltaT*0.5*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-2))/deltaR 
			+ 0.25*deltaT*deltaT*((-0.5*(3-gamma)*(bins[0]->U*bins[0]->U + bins[rgridNumber-1]->U*bins[rgridNumber-1]->U)*(densityFlux(0) - densityFlux(rgridNumber-1)) + 0.5*(3-gamma)*(bins[rgridNumber-1]->U*bins[rgridNumber-1]->U + bins[rgridNumber-2]->U*bins[rgridNumber-2]->U)*(densityFlux(rgridNumber-1) - densityFlux(rgridNumber-2))) 
			+ (2*(3-gamma)*(bins[0]->U + bins[rgridNumber-1]->U)*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1)) - 2*(3-gamma)*(bins[rgridNumber-1]->U+bins[rgridNumber-2]->U)*(fullMomentumFlux(rgridNumber-1) - fullMomentumFlux(rgridNumber-2))) 
			+ (2*(gamma - 1)*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-1)) - 2*(gamma - 1)*(fullEnergyFlux(rgridNumber-1) - fullEnergyFlux(rgridNumber-2))))/(deltaR*deltaR);
		newEnergy[rgridNumber-1] = newEnergy[rgridNumber-1] - deltaT*0.5*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-2))/deltaR 
			+ 0.25*deltaT*deltaT*(((-gamma*bins[0]->U*bins[0]->fullEnergy()/bins[0]->density + (gamma - 1)*bins[0]->U*bins[0]->U*bins[0]->U - gamma*bins[rgridNumber-1]->U*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density + (gamma - 1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U)*(densityFlux(0) - densityFlux(rgridNumber-1)) - (-gamma*bins[rgridNumber-1]->U*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density + (gamma - 1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U - gamma*bins[rgridNumber-2]->U*bins[rgridNumber-2]->fullEnergy()/bins[rgridNumber-2]->density + (gamma - 1)*bins[rgridNumber-2]->U*bins[rgridNumber-2]->U*bins[rgridNumber-2]->U)*(densityFlux(rgridNumber-1) - densityFlux(rgridNumber-2))) 
			+ ((gamma*bins[0]->fullEnergy()/bins[0]->density - 3*(gamma-1)*bins[0]->U*bins[0]->U/2 + gamma*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density - 3*(gamma-1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U/2)*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1)) - (gamma*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density - 3*(gamma-1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U/2 + gamma*bins[rgridNumber-2]->fullEnergy()/bins[rgridNumber-2]->density - 3*(gamma-1)*bins[rgridNumber-2]->U*bins[rgridNumber-2]->U/2)*(fullMomentumFlux(rgridNumber-1) - fullMomentumFlux(rgridNumber-2))) 
			+ (gamma*(bins[0]->U + bins[rgridNumber-1]->U)*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-1)) - gamma*(bins[rgridNumber-1]->U + bins[rgridNumber-2]->U)*(fullEnergyFlux(rgridNumber-1) - fullEnergyFlux(rgridNumber-2))))/(deltaR*deltaR);*/

	if(bins[rgridNumber - 1]->U >= 0){
		newDensity[rgridNumber-1] = newDensity[rgridNumber-1] - deltaT*(densityFlux(rgridNumber - 1) - densityFlux(rgridNumber-2))/deltaR 
			+ 0.25*deltaT*deltaT*((2*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1)) - 2*(fullMomentumFlux(rgridNumber-1) - fullMomentumFlux(rgridNumber-2))))/(deltaR*deltaR);
		newMomentum[rgridNumber-1] = newMomentum[rgridNumber-1] - deltaT*(fullMomentumFlux(rgridNumber - 1) - fullMomentumFlux(rgridNumber-2))/deltaR 
			+ 0.25*deltaT*deltaT*((-0.5*(3-gamma)*(bins[0]->U*bins[0]->U + bins[rgridNumber-1]->U*bins[rgridNumber-1]->U)*(densityFlux(0) - densityFlux(rgridNumber-1)) + 0.5*(3-gamma)*(bins[rgridNumber-1]->U*bins[rgridNumber-1]->U + bins[rgridNumber-2]->U*bins[rgridNumber-2]->U)*(densityFlux(rgridNumber-1) - densityFlux(rgridNumber-2))) 
			+ (2*(3-gamma)*(bins[0]->U + bins[rgridNumber-1]->U)*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1)) - 2*(3-gamma)*(bins[rgridNumber-1]->U+bins[rgridNumber-2]->U)*(fullMomentumFlux(rgridNumber-1) - fullMomentumFlux(rgridNumber-2))) 
			+ (2*(gamma - 1)*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-1)) - 2*(gamma - 1)*(fullEnergyFlux(rgridNumber-1) - fullEnergyFlux(rgridNumber-2))))/(deltaR*deltaR);
		newEnergy[rgridNumber-1] = newEnergy[rgridNumber-1] - deltaT*(fullEnergyFlux(rgridNumber - 1) - fullEnergyFlux(rgridNumber-2))/deltaR 
			+ 0.25*deltaT*deltaT*(((-gamma*bins[0]->U*bins[0]->fullEnergy()/bins[0]->density + (gamma - 1)*bins[0]->U*bins[0]->U*bins[0]->U - gamma*bins[rgridNumber-1]->U*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density + (gamma - 1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U)*(densityFlux(0) - densityFlux(rgridNumber-1)) - (-gamma*bins[rgridNumber-1]->U*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density + (gamma - 1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U - gamma*bins[rgridNumber-2]->U*bins[rgridNumber-2]->fullEnergy()/bins[rgridNumber-2]->density + (gamma - 1)*bins[rgridNumber-2]->U*bins[rgridNumber-2]->U*bins[rgridNumber-2]->U)*(densityFlux(rgridNumber-1) - densityFlux(rgridNumber-2))) 
			+ ((gamma*bins[0]->fullEnergy()/bins[0]->density - 3*(gamma-1)*bins[0]->U*bins[0]->U/2 + gamma*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density - 3*(gamma-1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U/2)*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1)) - (gamma*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density - 3*(gamma-1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U/2 + gamma*bins[rgridNumber-2]->fullEnergy()/bins[rgridNumber-2]->density - 3*(gamma-1)*bins[rgridNumber-2]->U*bins[rgridNumber-2]->U/2)*(fullMomentumFlux(rgridNumber-1) - fullMomentumFlux(rgridNumber-2))) 
			+ (gamma*(bins[0]->U + bins[rgridNumber-1]->U)*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-1)) - gamma*(bins[rgridNumber-1]->U + bins[rgridNumber-2]->U)*(fullEnergyFlux(rgridNumber-1) - fullEnergyFlux(rgridNumber-2))))/(deltaR*deltaR);
	} else {
		newDensity[rgridNumber-1] = newDensity[rgridNumber-1] - deltaT*(densityFlux(0) - densityFlux(rgridNumber - 1))/deltaR 
			+ 0.25*deltaT*deltaT*((2*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1)) - 2*(fullMomentumFlux(rgridNumber-1) - fullMomentumFlux(rgridNumber-2))))/(deltaR*deltaR);
		newMomentum[rgridNumber-1] = newMomentum[rgridNumber-1] - deltaT*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber - 1))/deltaR 
			+ 0.25*deltaT*deltaT*((-0.5*(3-gamma)*(bins[0]->U*bins[0]->U + bins[rgridNumber-1]->U*bins[rgridNumber-1]->U)*(densityFlux(0) - densityFlux(rgridNumber-1)) + 0.5*(3-gamma)*(bins[rgridNumber-1]->U*bins[rgridNumber-1]->U + bins[rgridNumber-2]->U*bins[rgridNumber-2]->U)*(densityFlux(rgridNumber-1) - densityFlux(rgridNumber-2))) 
			+ (2*(3-gamma)*(bins[0]->U + bins[rgridNumber-1]->U)*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1)) - 2*(3-gamma)*(bins[rgridNumber-1]->U+bins[rgridNumber-2]->U)*(fullMomentumFlux(rgridNumber-1) - fullMomentumFlux(rgridNumber-2))) 
			+ (2*(gamma - 1)*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-1)) - 2*(gamma - 1)*(fullEnergyFlux(rgridNumber-1) - fullEnergyFlux(rgridNumber-2))))/(deltaR*deltaR);
		newEnergy[rgridNumber-1] = newEnergy[rgridNumber-1] - deltaT*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber - 1))/deltaR 
			+ 0.25*deltaT*deltaT*(((-gamma*bins[0]->U*bins[0]->fullEnergy()/bins[0]->density + (gamma - 1)*bins[0]->U*bins[0]->U*bins[0]->U - gamma*bins[rgridNumber-1]->U*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density + (gamma - 1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U)*(densityFlux(0) - densityFlux(rgridNumber-1)) - (-gamma*bins[rgridNumber-1]->U*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density + (gamma - 1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U - gamma*bins[rgridNumber-2]->U*bins[rgridNumber-2]->fullEnergy()/bins[rgridNumber-2]->density + (gamma - 1)*bins[rgridNumber-2]->U*bins[rgridNumber-2]->U*bins[rgridNumber-2]->U)*(densityFlux(rgridNumber-1) - densityFlux(rgridNumber-2))) 
			+ ((gamma*bins[0]->fullEnergy()/bins[0]->density - 3*(gamma-1)*bins[0]->U*bins[0]->U/2 + gamma*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density - 3*(gamma-1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U/2)*(fullMomentumFlux(0) - fullMomentumFlux(rgridNumber-1)) - (gamma*bins[rgridNumber-1]->fullEnergy()/bins[rgridNumber-1]->density - 3*(gamma-1)*bins[rgridNumber-1]->U*bins[rgridNumber-1]->U/2 + gamma*bins[rgridNumber-2]->fullEnergy()/bins[rgridNumber-2]->density - 3*(gamma-1)*bins[rgridNumber-2]->U*bins[rgridNumber-2]->U/2)*(fullMomentumFlux(rgridNumber-1) - fullMomentumFlux(rgridNumber-2))) 
			+ (gamma*(bins[0]->U + bins[rgridNumber-1]->U)*(fullEnergyFlux(0) - fullEnergyFlux(rgridNumber-1)) - gamma*(bins[rgridNumber-1]->U + bins[rgridNumber-2]->U)*(fullEnergyFlux(rgridNumber-1) - fullEnergyFlux(rgridNumber-2))))/(deltaR*deltaR);
	}
}


