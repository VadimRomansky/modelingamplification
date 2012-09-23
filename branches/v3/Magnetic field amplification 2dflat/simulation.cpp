#include "stdafx.h"
#include "simulation.h"
#include "SpaceBin.h"
#include "particle.h"
#include "output.h"
#include "util.h"
#include <omp.h>
//#include <Windows.h>

Simulation::Simulation(){
	allParticlesNumber = 0;
	energy = 0;
	theorEnergy = 0;
	momentumX = 0;
	theorMomentumX = 0;
	A = 1;
	Z = 1;
	startPDF = std::list <Particle*>();
	timeStep = defaultTimeStep;
}

Simulation::~Simulation(){
}

void Simulation::initializeProfile(){
	deltaX = (X2 - X1)/(xgridNumber );
	deltaY = (Y2 - Y1)/(ygridNumber );

	double X = X1 + deltaX/2;
	double Y = Y1 + deltaY/2;
	/*if(freeTimeEvaluatorType == 1){
		freeTimeEvaluator = new ConstantFreeTimeEvaluator(B0);
	} else if(freeTimeEvaluatorType == 2) {
		freeTimeEvaluator = new LinearFreeTimeEvaluator(B0,particleLocalMomentum);
	} else {
		freeTimeEvaluator = NULL;
	}*/
	bins = new SpaceBin**[xgridNumber];
    for(int i = 0; i < xgridNumber; ++i){
		Y = Y1 + deltaY/2;
		bins[i] = new SpaceBin*[ygridNumber];
		for(int j = 0; j < ygridNumber; ++j){
			double density = density0;
			double u;
			if((i+0.5)*(i+0.5) + (j+0.5)*(j+0.5) < (shockWaveIndex+0.5)*(shockWaveIndex+0.5)){
				u = U0;
			} else {
				u = 0;
			}
			double ux = u*(i + 0.5)/(sqrt(1.0*((i+0.5)*(i+0.5) + (j+0.5)*(j+0.5))));
			double uy = u*(j + 0.5)/(sqrt(1.0*((i+0.5)*(i+0.5) + (j+0.5)*(j+0.5))));
			bins[i][j] = new SpaceBin(X, deltaX, Y, deltaY, ux, uy, density, temperature, B0, i, j, smallAngleScattering, freeTimeEvaluationType);
			Y = Y + deltaY;
		}
		X = X + deltaX;
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

				int* index = SpaceBin::binByCoordinates(particle->absoluteX, particle->absoluteY, X1, deltaX, xgridNumber, Y1, deltaY, ygridNumber);

				double time = 0;
				if((index[0] == -1) || (index[1] == -1)){
					printf("(index[0] == -1) || (index[1] == -1)\n");
					printf("10000000*x = %lf 100000000*y = %lf\n", 100000000*particle->absoluteX, 10000000*particle->absoluteY);
				}
				while((index[0] < xgridNumber) && (index[1] < ygridNumber) && (index[0] >= 0) && (index[1] >= 0)){
					j++;
					SpaceBin* bin = bins[index[0]][index[1]];
					if(time != time){
						printf("time != time\n");
					}
					int* tempIndex = bin->propagateParticle(particle,time, timeStep, xgridNumber, ygridNumber);
					delete[] index;
					index = tempIndex;
					if((index[0] < -1) || (index[1] < -1)){
						printf("(index[0] < -1) || (index[1] < -1)\n");
					}
					if((index[0] < 0) || (index[1] < 0)){
						if(index[0] == -1){
							index[0] = 0;
							particle->absoluteMomentumX= - particle->absoluteMomentumX;
							particle->moveToBinLeftX(bin);
							theorMomentumX += 2*particle->absoluteMomentumX*particle->weight;
						}
						if(index[1] == -1){
							index[1] = 0;
							particle->absoluteMomentumY= - particle->absoluteMomentumY;
							particle->moveToBinLeftY(bin);
							theorMomentumY += 2*particle->absoluteMomentumY*particle->weight;
						}
					}
					if (time >= timeStep){
						break;
					}
				}
				delete[] index;
			}

			//printf("%s", "collecting velocity, density and crflux\n");
			collectAverageVelocity();
			//printf("%s", "removing escaped particles\n");
			removeEscapedParticles();
			//printf("%s","updating energy\n");
			//updateCosmicRayBoundMomentum(itNumber % writeParameter == 0);

			//outputEnergyPDF(introducedParticles,"./output/tamc_energy_pdf.dat");
			resetDetectors();
			//printf("%s","iteration  ");
			//printf("%d\n",itNumber);
		}
		updateEnergy();
		if(itNumber % writeParameter == 0){
			outputProfile(bins, xgridNumber, ygridNumber);
			outputEnergyPDF(introducedParticles,"./output/tamc_energy_pdf.dat");
			//findShockWavePoint();
			printf("%s","iteration  ");
			printf("%d\n",itNumber);
			printf("outputing\n");
			outputParticles(introducedParticles,"./output/particles.dat");
			outputPDF(introducedParticles,"./output/tamc_pdf.dat");
			if(escapedParticles.size() > 10){
				outputPDF(escapedParticles,"./output/escaped_pdf.dat");
			}
			outputEnergyPDF(introducedParticles,"./output/tamc_energy_pdf.dat");
			outIteration = fopen("./output/tamc_iteration.dat","a");
			radialFile = fopen("./output/tamc_radial_profile.dat","a");
			fprintf(outIteration,"%d %lf %lf %lf %lf %d %lf %lf\n",itNumber, energy, theorEnergy, momentumX, theorMomentumX, introducedParticles.size(), particlesWeight, shockWavePoint);
			fclose(outIteration);
			outputRadialProfile(bins[0],radialFile, ygridNumber);
			fclose(radialFile);
		}
	}
}

////   
void Simulation::resetDetectors(){
	for(int i = 0; i < xgridNumber; ++i){
		for(int j = 0; j < ygridNumber; ++j){
			SpaceBin* bin = bins[i][j];
			bin->resetDetectors();			
		}	
	}
}

std::vector <Particle*> Simulation::getParticles(){
	std::vector<Particle*> list = std::vector<Particle*>();
	int cosmicRayNumber = 0;
	int particleBinNumber = 0;
	allParticlesNumber = 0;
	for(int i = 0; i < xgridNumber; ++i){
		for(int j = 0; j < ygridNumber; ++j){
			particleBinNumber = 0;
			SpaceBin* bin = bins[i][j];
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
				momentumY += particle->absoluteMomentumY*particle->weight;

			}
		}
	}
	theorEnergy = energy;
	theorMomentumX = momentumX;
	theorMomentumY = momentumY;
	return list;
}

void Simulation::removeEscapedParticles(){
	std::vector<Particle*> list = std::vector<Particle*>();
	std::vector<Particle*>::iterator it = introducedParticles.begin();
	while(it != introducedParticles.end()){
		Particle* particle = *it;
		++it;
		if((particle->absoluteX < X1) || (particle->absoluteX > X2) || (particle->absoluteY < Y1) || (particle->absoluteY > Y2)){
			theorEnergy -= particle->getEnergy()*particle->weight;
			theorMomentumX -= particle->absoluteMomentumX*particle->weight;
			theorMomentumY -= particle->absoluteMomentumY*particle->weight;
			escapedParticles.push_back(particle);
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
	double** count = new double*[xgridNumber];
	double** averageVelocityX = new double*[xgridNumber];
	double** averageVelocityY = new double*[xgridNumber];
	double** absoluteEnergy = new double*[xgridNumber];
	double** absoluteTheorEnergy = new double*[xgridNumber];

	for(int i = 0; i < xgridNumber; ++i){
		count[i] = new double[ygridNumber];
		averageVelocityX[i] = new double[ygridNumber];
		averageVelocityY[i] = new double[ygridNumber];
		absoluteEnergy[i] = new double[ygridNumber];
		absoluteTheorEnergy[i] = new double[ygridNumber];
		for(int j = 0; j < ygridNumber; ++j){
			count[i][j] = 0;
			averageVelocityX[i][j] = 0;
			averageVelocityY[i][j] = 0;
			absoluteEnergy[i][j] = 0;
			absoluteTheorEnergy[i][j] = 0;
			bins[i][j]->density = 0;
		}
	}

	std::vector<Particle*>::iterator it = introducedParticles.begin();

	while(it != introducedParticles.end()){
		Particle* particle = *it;
		++it;
		double theta;
		double r = particle->absoluteX;

		int* index = SpaceBin::binByCoordinates(particle->absoluteX, particle->absoluteY, X1, deltaX, xgridNumber, Y1, deltaY, ygridNumber);
		if((index[0] >= 0)&&(index[0] < xgridNumber) && (index[1] >= 0)&&(index[1] < xgridNumber)){
			count[index[0]][index[1]] += particle->weight;
			bins[index[0]][index[1]]->density += particle->mass*particle->weight;
			bins[index[0]][index[1]]->particleMomentaZ.push_back(particle->localMomentumX);
			bins[index[0]][index[1]]->particleWeights.push_back(particle->weight);
		}
	}

	/*it = introducedParticles.begin();

	while(it != introducedParticles.end()){
		Particle* particle = *it;
		++it;
		int* index = SpaceBin::binByCoordinates(particle->absoluteX, particle->absoluteY, X1, deltaX, xgridNumber, Y1, deltaY, ygridNumber);
		if((index[0] >= 0)&&(index[0] < xgridNumber) && (index[1] >= 0) && (index[1] < ygridNumber)){
			double energy = particle->getEnergy();
			double deltaP;
			if(count[index[0]][index[1]] > epsilon){
				deltaP = bins[index[0]][index[1]]->momentumDifference/count[index[0]][index[1]];
			} else {
				deltaP = 0;
			}
			if(particle->absoluteMomentum*particle->absoluteMomentum <= particle->absoluteMomentumX*particle->absoluteMomentumX){
				particle->absoluteMomentumX -= deltaP;
				particle->absoluteMomentum = particle->absoluteMomentumX;
			} else {
				double particleMomentumY = sqrt(particle->absoluteMomentum*particle->absoluteMomentum - particle->absoluteMomentumX*particle->absoluteMomentumX);
				particle->absoluteMomentumX -= deltaP;
				particle->absoluteMomentum = sqrt(particle->absoluteMomentumX*particle->absoluteMomentumX + particleMomentumY*particleMomentumY);
			}
			bins[index[0]][index[1]]->energyDifference += (particle->getEnergy() - energy)*particle->weight;
			absoluteEnergy[index[0]][index[1]] +=particle->getEnergy()*particle->weight;
		}
	}*/

	it = introducedParticles.begin();

	while(it != introducedParticles.end()){
		Particle* particle = *it;
		++it;
		double theta;
		double r = particle->absoluteX;

		int* index = SpaceBin::binByCoordinates(particle->absoluteX, particle->absoluteY, X1, deltaX, xgridNumber, Y1, deltaY, ygridNumber);
		if((index[0] >= 0)&&(index[0] < xgridNumber) && (index[1] >= 0) && (index[1] < ygridNumber)){
			averageVelocityX[index[0]][index[1]] += particle->getAbsoluteVX()*particle->weight;
			averageVelocityY[index[0]][index[1]] += particle->getAbsoluteVY()*particle->weight;
		}
	}

	/*for(int i = 0; i < xgridNumber; ++i){
		absoluteTheorEnergy[i] = absoluteEnergy[i] - bins[i]->energyDifference;
		if(absoluteTheorEnergy[i] <= 0){
			absoluteTheorEnergy[i] = absoluteEnergy[i];
		}
	}*/

	for(int i = 0; i < xgridNumber; ++i){
		for(int j = 0; j < ygridNumber; ++j){
			bins[i][j]->density /= bins[i][j]->volume;
			if(count[i][j] > epsilon){
				averageVelocityX[i][j] /= count[i][j];
				averageVelocityY[i][j] /= count[i][j];
				bins[i][j]->averageVelocityX = averageVelocityX[i][j];
				bins[i][j]->averageVelocityY = averageVelocityY[i][j];
				bins[i][j]->Ux = averageVelocityX[i][j];
				bins[i][j]->Uy = averageVelocityY[i][j];
				//bins[i]->U = averageVelocity[i] - bins[i]->momentumDifference/(massProton*count[i]);
			} else {
				printf("0 particles in bin\n");
				bins[i][j]->averageVelocityX = 0;
				bins[i][j]->averageVelocityY = 0;
				bins[i][j]->Ux = 0;
				bins[i][j]->Uy = 0;
			}
		}
	}

	/*it = introducedParticles.begin();

	while(it != introducedParticles.end()){
		Particle* particle = *it;
		++it;
		int index = SpaceBin::binByCoordinates(particle->absoluteX,upstreamR,deltaR, rgridNumber);
		if((index >= 0)&&(index < rgridNumber)){
			double localEnergy = absoluteEnergy[index] - bins[index]->density*bins[index]->volume*(bins[index]->U*bins[index]->U)/2;
			double localTheorEnergy = absoluteTheorEnergy[index] - bins[index]->density*bins[index]->volume*(bins[index]->U*bins[index]->U)/2;
			if(localEnergy < 0){
				printf("localEnergy < 0\n");
			}
			if(localTheorEnergy < 0){
				printf("localTheorEnergy < 0\n");
			}
			double momentumRelation = sqrt(localTheorEnergy/localEnergy);
			particle->localMomentum *= momentumRelation;
			particle->localMomentumX *= momentumRelation;
			particle->setAbsoluteMomentum(bins[index]->U);
		}
	}*/

	for(int i = 0; i < xgridNumber; ++i){
		delete[] count[i];
		delete[] averageVelocityX[i];
		delete[] averageVelocityY[i];
		delete[] absoluteEnergy[i];
		delete[] absoluteTheorEnergy[i];
	}

	delete[] count;
	delete[] averageVelocityX;
	delete[] averageVelocityY;
	delete[] absoluteEnergy;
	delete[] absoluteTheorEnergy;
}

void Simulation::updateCosmicRayBoundMomentum(bool write){
	for(int i = 0; i < xgridNumber; ++i){
		for(int j = 0; j < ygridNumber; ++j){
			bins[i][j]->updateCosmicRayBoundMomentum(write);
		}
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
		momentumY += particle->absoluteMomentumY*particle->weight;
		particlesWeight += particle->weight;
		++it;
	}
}

/*void Simulation::findShockWavePoint(){
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
}*/