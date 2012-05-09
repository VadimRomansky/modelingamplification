#include "stdafx.h"
#include <list>
#include "output.h"
#include "simulation.h"
#include "constants.h"
#include "particle.h"
#include "util.h"

void output(Simulation& simulation){
	double maxp, minp;
	/*FILE* outProfile = fopen("./output/tamc_grid.dat","w");
	//FILE* outPDF = fopen("tamc_pdf.dat","w");
	double massFlux0 = simulation.xbins[0]->massFlux;
	double momentaFlux0 = simulation.xbins[0]->momentaFlux;
	double energyFlux0 = simulation.xbins[0]->energyFlux;
	for(int i = 0; i < xgridNumber; ++i){
		Xbin* bin = simulation.xbins[i];
		fprintf(outProfile, "%lf %lf %lf %lf %1.5lf %1.5lf %1.5lf",(bin->left+bin->right)/2, bin->U, bin->density, bin->temperature, bin->massFlux/massFlux0, bin->momentaFlux/momentaFlux0, bin->energyFlux/energyFlux0);
		fprintf(outProfile,"%s","\n");
	}
	fclose(outProfile);*/
    SpaceBin* bin = simulation.bins[0][0][0];
	simulation.updateMaxMinP(maxp,minp);
	outputPDF(simulation.startPDF,"./output/tamc_pdf_start.dat");
	outputPDF(bin->detectedParticlesR2,"./output/tamc_pdf_FEB.dat");
	//outputTurbulenceSpectrum(bin->magneticField,"./output/tamc_turb_FEB.dat",simulation.minK,simulation.maxK);
	bin = simulation.bins[simulation.shockWavePoint][0][0];
	outputPDF(bin->detectedParticlesR2,"./output/tamc_pdf_zero.dat");
	//outputTurbulenceSpectrum(bin->magneticField,"./output/tamc_turb_zero.dat",simulation.minK,simulation.maxK);
	bin = simulation.bins[rgridNumber-1][0][0];
	outputPDF(bin->detectedParticlesR2,"./output/tamc_pdf_down.dat");
	//outputTurbulenceSpectrum(bin->magneticField,"./output/tamc_turb_down.dat",simulation.minK,simulation.maxK);
	//outputMagneticField(simulation.xbins,"./output/tamc_field.dat");
	//outputParticlePath(simulation.xbins[xgridNumber - 1]->detectedParticlesRight,"./output/tamc_cosmic_ray_path.dat","./output/tamc_not_cosmic_ray_path.dat");
	//outputParticlesCount(simulation.xbins,"./output/tamc_particle_count.dat");
	//FILE* outParticles = fopen("./output/tamc_particle_count.dat","w");
	//fprintf(outParticles,"%d \n",simulation.allParticlesNumber);
	/*for(int i = 0; i < xgridNumber; ++i){
		fprintf(outParticles,"%lf %d\n", (simulation.xbins[i]->left+simulation.xbins[i]->right)/2,simulation.xbins[i]->detectedParticlesLeft.size());
		fprintf(outParticles,"%lf %d\n", (simulation.xbins[i]->left+simulation.xbins[i]->right)/2,simulation.xbins[i]->detectedParticlesRight.size());
	}*/
	//fclose(outParticles);
}

void outputPDF(std::list< Particle*>& list,const char* fileName){
	FILE* outPDF = fopen(fileName,"w");
	std::list<Particle*>::iterator it = list.begin();
	double minp;
	double maxp;
	minp = (*it)->absoluteMomentum;
	maxp = minp;
	++it;
	while(it != list.end()){
		Particle* particle = *it;
		if(maxp < particle->absoluteMomentum){
			maxp = particle->absoluteMomentum;
		}
		if(maxp < particle->initialMomentum){
			maxp = particle->initialMomentum;
		}
		if(minp > particle->absoluteMomentum){
			minp = particle->absoluteMomentum;
		}
		if(minp > particle->initialMomentum){
			minp = particle->initialMomentum;
		}
		++it;
	}
	double* distribution = new double[pgridNumber];
	double* startDistribution = new double[pgridNumber];
	for (int i = 0; i < pgridNumber; ++i){
		distribution[i] = 0;
		startDistribution[i] = 0;
	}
	int particleNumber = 0;
	double mass;
	double deltap = (maxp - minp)/(pgridNumber - 1);
	it = list.begin();
	while(it != list.end()){
		particleNumber++;
		Particle* particle = *it;
		mass = particle->mass;
		++it;
		double p = particle->absoluteMomentum;
		if (p >= maxp){
			distribution[pgridNumber-1]+=particle->weight;
		}
		for(int i =0; i< pgridNumber-1; ++i){
			if (p <minp + deltap*(i + 1)){
				distribution[i] += particle->weight;
				break;
			}
		}
		p = particle->initialMomentum;
		if (p >= maxp){
			startDistribution[pgridNumber-1]+=particle->weight;
		}
		for(int i =0; i< pgridNumber-1; ++i){
			if (p <minp + deltap*(i + 1)){
				startDistribution[i] += particle->weight;
				break;
			}
		}
	}
	for(int i = 0; i < pgridNumber; ++i){
		fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(minp + i*deltap), (distribution[i]), (startDistribution[i]));
	}
	delete[] distribution;
	delete[] startDistribution;
	fclose(outPDF); 
}

void outputPDF(std::vector< Particle*>& list,const char* fileName){
	FILE* outPDF = fopen(fileName,"w");
	std::vector<Particle*>::iterator it = list.begin();
	double minp;
	double maxp;
	minp = (*it)->absoluteMomentum;
	maxp = minp;
	++it;
	while(it != list.end()){
		Particle* particle = *it;
		if(maxp < particle->absoluteMomentum){
			maxp = particle->absoluteMomentum;
		}
		if(maxp < particle->initialMomentum){
			maxp = particle->initialMomentum;
		}
		if(minp > particle->absoluteMomentum){
			minp = particle->absoluteMomentum;
		}
		if(minp > particle->initialMomentum){
			minp = particle->initialMomentum;
		}
		++it;
	}
	double* distribution = new double[pgridNumber];
	double* startDistribution = new double[pgridNumber];
	for (int i = 0; i < pgridNumber; ++i){
		distribution[i] = 0;
		startDistribution[i] = 0;
	}
	int particleNumber = 0;
	double mass;
	double deltap = (maxp - minp)/(pgridNumber - 1);
	it = list.begin();
	while(it != list.end()){
		particleNumber++;
		Particle* particle = *it;
		mass = particle->mass;
		++it;
		double p = particle->absoluteMomentum;
		if (p >= maxp){
			distribution[pgridNumber-1]+=particle->weight;
		}
		for(int i =0; i< pgridNumber-1; ++i){
			if (p <minp + deltap*(i + 1)){
				distribution[i] += particle->weight;
				break;
			}
		}
		p = particle->initialMomentum;
		if (p >= maxp){
			startDistribution[pgridNumber-1]+=particle->weight;
		}
		for(int i =0; i< pgridNumber-1; ++i){
			if (p <minp + deltap*(i + 1)){
				startDistribution[i] += particle->weight;
				break;
			}
		}
	}
	for(int i = 0; i < pgridNumber; ++i){
		fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(minp + i*deltap), (distribution[i]), (startDistribution[i]));
	}
	delete[] distribution;
	delete[] startDistribution;
	fclose(outPDF); 
}

void outputEnergyPDF(std::vector< Particle*>& list,const char* fileName){
	if(list.size() > 0){
		FILE* outPDF = fopen(fileName,"w");
		std::vector<Particle*>::iterator it = list.begin();
		double mine;
		double maxe;
		mine = (*it)->getEnergy();
		maxe = mine;
		++it;
		while(it != list.end()){
			Particle* particle = *it;
			if(maxe < particle->getEnergy()){
				maxe = particle->getEnergy();
			}
			if(maxe < particle->getInitialEnergy()){
				maxe = particle->getInitialEnergy();
			}
			if(mine > particle->getEnergy()){
				mine = particle->getEnergy();
			}
			if(mine > particle->getInitialEnergy()){
				mine = particle->getInitialEnergy();
			}
			++it;
		}
		if(mine < 0){
			mine = epsilon*maxe;
		}
		double* distribution = new double[pgridNumber];
		double* startDistribution = new double[pgridNumber];
		for (int i = 0; i < pgridNumber; ++i){
			distribution[i] = 0;
			startDistribution[i] = 0;
		}
		int particleNumber = 0;
		double mass;
		double deltae = (log(maxe) - log(mine))/(pgridNumber - 1);
		it = list.begin();
		while(it != list.end()){
			particleNumber++;
			Particle* particle = *it;
			mass = particle->mass;
			++it;
			double p = particle->getEnergy();
			if (log(p) >= log(maxe)){
				distribution[pgridNumber-1]+=particle->weight;
			}
			for(int i =0; i< pgridNumber-1; ++i){
				if (log(p) <log(mine) + deltae*(i + 1)){
					distribution[i] += particle->weight;
					break;
				}
			}
			p = particle->getInitialEnergy();
			if (log(p) >= log(maxe)){
				startDistribution[pgridNumber-1]+=particle->weight;
			}
			for(int i =0; i< pgridNumber-1; ++i){
				if (log(p) <log(mine) + deltae*(i + 1)){
					startDistribution[i] += particle->weight;
					break;
				}
			}
		}
		for(int i = 0; i < pgridNumber; ++i){
			fprintf(outPDF,"%lf %lf %lf\n", (log(mine) + i*deltae)/log(10.0), (distribution[i])/exp(log(mine) + i*deltae), (startDistribution[i])/exp(log(mine) + i*deltae));
		}
		delete[] distribution;
		delete[] startDistribution;
		fclose(outPDF); 
	}
}

void outputEnergyPDF(std::list< Particle*>& list,const char* fileName){
	if(list.size() > 0){
		FILE* outPDF = fopen(fileName,"w");
		std::list<Particle*>::iterator it = list.begin();
		double mine;
		double maxe;
		mine = (*it)->getEnergy();
		maxe = mine;
		++it;
		while(it != list.end()){
			Particle* particle = *it;
			if(maxe < particle->getEnergy()){
				maxe = particle->getEnergy();
			}
			if(maxe < particle->getInitialEnergy()){
				maxe = particle->getInitialEnergy();
			}
			if(mine > particle->getEnergy()){
				mine = particle->getEnergy();
			}
			if(mine > particle->getInitialEnergy()){
				mine = particle->getInitialEnergy();
			}
			++it;
		}
		if(mine < 0){
			mine = epsilon*maxe;
		}
		double* distribution = new double[pgridNumber];
		double* startDistribution = new double[pgridNumber];
		for (int i = 0; i < pgridNumber; ++i){
			distribution[i] = 0;
			startDistribution[i] = 0;
		}
		int particleNumber = 0;
		double mass;
		double deltae = (log(maxe) - log(mine))/(pgridNumber - 1);
		it = list.begin();
		while(it != list.end()){
			particleNumber++;
			Particle* particle = *it;
			mass = particle->mass;
			++it;
			double p = particle->getEnergy();
			if (log(p) >= log(maxe)){
				distribution[pgridNumber-1]+=particle->weight;
			}
			for(int i =0; i< pgridNumber-1; ++i){
				if (log(p) <log(mine) + deltae*(i + 1)){
					distribution[i] += particle->weight;
					break;
				}
			}
			p = particle->getInitialEnergy();
			if (log(p) >= log(maxe)){
				startDistribution[pgridNumber-1]+=particle->weight;
			}
			for(int i =0; i< pgridNumber-1; ++i){
				if (log(p) <log(mine) + deltae*(i + 1)){
					startDistribution[i] += particle->weight;
					break;
				}
			}
		}
		for(int i = 0; i < pgridNumber; ++i){
			fprintf(outPDF,"%lf %lf %lf\n", (log(mine) + i*deltae)/log(10.0), (distribution[i])/exp(log(mine) + i*deltae), (startDistribution[i])/exp(log(mine) + i*deltae));
		}
		delete[] distribution;
		delete[] startDistribution;
		fclose(outPDF); 
	}
}

void outputRPDF(std::list< Particle*>& list,const char* fileName){
	FILE* outPDF = fopen(fileName,"w");
	std::list<Particle*>::iterator it = list.begin();
	double maxp;
	maxp = (*it)->localMomentumZ;
	++it;
	double weight = 0;
	weight += (*it)->weight;
	while(it != list.end()){
		Particle* particle = *it;
		weight += particle->weight;
		if(maxp < abs(particle->localMomentumZ)){
			maxp = abs(particle->localMomentumZ);
		}
		++it;
	}
	double* distribution = new double[pgridNumber];
	double* startDistribution = new double[pgridNumber];
	for (int i = 0; i < pgridNumber; ++i){
		distribution[i] = 0;
		startDistribution[i] = 0;
	}
	int particleNumber = 0;
	double mass;
	double deltap = (2*maxp)/(pgridNumber);
	it = list.begin();
	while(it != list.end()){
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
	double temperature;
	double sigma;
	int zeroP = pgridNumber/2;
	int sigmaIndex;
	double maxDistribution = (distribution[zeroP] + distribution[zeroP + 1])/2;
	for(int i = zeroP;i < pgridNumber; ++i){
		if(distribution[i] < maxDistribution/exp(0.5)){
			sigmaIndex = i - 1;
			break;
		}
	}
	double d1 = distribution[sigmaIndex];
	double d2 = distribution[sigmaIndex + 1];

	sigma =  (maxDistribution/exp(0.5) - d1)*deltap/(d2 - d1) + (sigmaIndex + 0.5 - zeroP)*deltap;
	temperature = sigma*sigma/(massProton*kBoltzman);
	for(int i = 0; i < pgridNumber; ++i){
		fprintf(outPDF,"%lf %lf %lf\n", 10000000000000000000.0*(-maxp + i*deltap), distribution[i], maxwell(-maxp + (i + 1/2)*deltap, massProton, temperature));
	}
	delete[] distribution;
	delete[] startDistribution;
	fclose(outPDF); 
}

void outputStartPDF(std::list< Particle*>& l,const char* fileName, Simulation& simulation,double minp,double maxp){
	FILE* outPDF = fopen(fileName,"w");
	std::list<Particle*>::iterator it = l.begin();
	double deltap = (maxp - minp)/(pgridNumber - 1);
	double* distribution = new double[pgridNumber];
	for (int i = 0; i < pgridNumber; ++i){
		distribution[i] = 0;
	}
	int particleNumber = 0;
	double mass;
	while(it != l.end()){
		particleNumber++;
		Particle* particle = *it;
		mass = particle->mass;
		++it;
		//double p = particle.momentum - particle.mass*simulation.U0;
		double p = particle->initialMomentum;
		//if (p >= maxMomentum){
		if (p >= maxp){
			distribution[pgridNumber-1]+=particle->weight;
		}
		for(int i =0; i< pgridNumber-1; ++i){
			//if (p < maxMomentum*(i+1)/(pgridNumber-2))	{
			if (p <minp + deltap*(i + 1)){
				distribution[i] += particle->weight;
				break;
			}
		}
	}
	for(int i = 0; i < pgridNumber; ++i){
		fprintf(outPDF,"%lf %lf\n", 100000000000*(minp + i*deltap), (distribution[i]));
	}
	delete[] distribution;
	fclose(outPDF); 
}

void outputTurbulenceSpectrum(double* w, const char* fileName, double minK, double maxK){
	FILE* outTurb = fopen(fileName,"w");
	for(int i = 0; i < kgridNumber; ++i){
		fprintf(outTurb,"%lf %lf\n", 100000000000000*(minK + i*(maxK - minK)/(kgridNumber-1)),100000000000000*w[i]);
	}
	fclose(outTurb);
}

void outputMagneticField(SpaceBin**** bins, const char* fileName){
	FILE* outField = fopen(fileName,"w");
	for(int i = 0; i < rgridNumber; ++i){
		fprintf(outField,"%lf %lf\n", (bins[i][0][0]->r2+bins[i][0][0]->r1)/2,bins[i][0][0]->B);
	}
	fclose(outField);
}

void outputParticlePath(std::list<Particle*>& list,const char* cosmicRayFileName,const char* notCosmicRayFileName){
	/*FILE* file1 = fopen(cosmicRayFileName,"w");
	FILE* file2 = fopen(notCosmicRayFileName,"w");
	std::list<Particle*>::iterator it = list.begin();
	Particle cosmicRayParticle;
	Particle notCosmicRayParticle;
	bool findCosmicRay = false;
	bool findNotCosmicRay = false;
	while(((!findCosmicRay) || (!findNotCosmicRay)) && (it != list.end())){
		Particle* particle = *it;
		if ((particle->isCosmicRay) && (!findCosmicRay)){
			cosmicRayParticle = *particle;
			findCosmicRay = true;
		}

		if ((!particle->isCosmicRay) && (!findNotCosmicRay)){
			notCosmicRayParticle = *particle;
			findNotCosmicRay = true;
		}

		if ((findCosmicRay) && (findNotCosmicRay)){
			break;
		}

		++it;

	}

	std::list<double>::iterator pathIterator = cosmicRayParticle.path.begin();
	fprintf(file1,"%d %s",cosmicRayParticle.path.size(),"\n");
	while (pathIterator != cosmicRayParticle.path.end()){
		fprintf(file1,"%lf %s",*pathIterator,"\n");
		++pathIterator;
	}
	fclose(file1);

	pathIterator = notCosmicRayParticle.path.begin();
	fprintf(file2,"%d %s",notCosmicRayParticle.path.size(),"\n");
	while (pathIterator != notCosmicRayParticle.path.end()){
		fprintf(file2,"%lf %s",*pathIterator,"\n");
		++pathIterator;
	}
	fclose(file2);*/
}

void outputRadialProfile(SpaceBin**** bins, int thetaNumber, int phiNumber, FILE* outProfile){
	//FILE* outPDF = fopen("tamc_pdf.dat","w");
	//double massFlux0 = simulation.xbins[0]->massFlux;
	//double momentaFlux0 = simulation.xbins[0]->momentaFlux;
	//double energyFlux0 = simulation.xbins[0]->energyFlux;
	for(int i = 0; i < rgridNumber; ++i){
		SpaceBin* bin = bins[i][thetaNumber][phiNumber];
		fprintf(outProfile, "%lf %lf %lf %lf %lf %lf",bin->r, bin->U, bin->averageVelocity, 100000*bin->density, bin->temperature);
		fprintf(outProfile,"%s","\n");
	}
	//fclose(outProfile);
}

void outputThetaProfile(SpaceBin**** bins, int rNumber, int phiNumber){
	FILE* outProfile = fopen("./output/tamc_grid_theta.dat","w");
	//FILE* outPDF = fopen("tamc_pdf.dat","w");
	//double massFlux0 = simulation.xbins[0]->massFlux;
	//double momentaFlux0 = simulation.xbins[0]->momentaFlux;
	//double energyFlux0 = simulation.xbins[0]->energyFlux;
	for(int i = 0; i < thetagridNumber; ++i){
		SpaceBin* bin = bins[rNumber][i][phiNumber];
		fprintf(outProfile, "%lf %lf %lf %lf",bin->theta, bin->U, bin->density, bin->temperature);
		fprintf(outProfile,"%s","\n");
	}
	fclose(outProfile);
}

void outputPhiProfile(SpaceBin**** bins, int rNumber, int thetaNumber){
	FILE* outProfile = fopen("./output/tamc_grid_phi.dat","w");
	//FILE* outPDF = fopen("tamc_pdf.dat","w");
	//double massFlux0 = simulation.xbins[0]->massFlux;
	//double momentaFlux0 = simulation.xbins[0]->momentaFlux;
	//double energyFlux0 = simulation.xbins[0]->energyFlux;
	for(int i = 0; i < phigridNumber; ++i){
		SpaceBin* bin = bins[rNumber][thetaNumber][i];
		fprintf(outProfile, "%lf %lf %lf %lf",bin->phi, bin->U, bin->density, bin->temperature);
		fprintf(outProfile,"%s","\n");
	}
	fclose(outProfile);
}

void outputShockWave(std::list<double> points, std::list<double> velocity){
	FILE* outWave = fopen("./output/tamc_shock_wave.dat", "w");
	std::list<double>::iterator it1 = points.begin();
	std::list<double>::iterator it2 = velocity.begin();
	int i = 0;
	while( (it1 != points.end()) && (it2 != velocity.end())){
		double velocity = *it2;
		double point = *it1;
		fprintf(outWave,"%d %lf %lf \n",i, point, velocity); 
		++it1;
		++it2;
		++i;
	}
	fclose(outWave);
}

void outputAverageVz(double minp, double maxp, double* averageVz, const char* fileName){
	FILE* outVz = fopen(fileName, "w");
	double deltap = (maxp - minp)/pgridNumber;
	for(int i = 0; i < pgridNumber; ++i){
		fprintf(outVz,"%lf %lf \n",1000000000000000000000000.0*(minp + i *deltap), averageVz[i]);
	}
	fclose(outVz);
}

void outputParticles(std::vector<Particle*>& particles, const char* fileName){
	FILE* outParticles = fopen(fileName, "w");
	std::vector<Particle*>::iterator it = particles.begin();
	while(it != particles.end()){
		fprintf(outParticles ,"%lf %lf \n",(*it)->getAbsoluteR(), (*it)->absoluteMomentum*1000000000000000000000000.0);
		++it;
	}
	fclose(outParticles );
}

/*void outputPartclesCount(Xbin** bins, const char* fileName){
	FILE* outParticles = fopen(fileName,"w");
	for(int i = 0; i < xgridNumber; ++i){
		fprintf(outParticles,"%lf %d\n", (bins[i]->left+bins[i]->right)/2,bins[i]->detectedParticles.size);
	}
	fclose(outParticles);
}*/