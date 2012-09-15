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
	//sortedParticles = new std::list<Particle*>[kgridNumber];
}
SpaceBin::SpaceBin(double R, double deltar, double u, double rho, double t, double b, int i, bool scattering){
	r = R;
	r1 = r - deltar/2;
	r2 = r + deltar/2;

	volume = deltar;

	U = u;

	averageVelocity = U;

	density = rho;
	temperature = t;

	B0 = b;

	numberR = i;

	smallAngleScattering = scattering;


	particles.clear();

	//TODO!!!!!!!!!!!!!!!!
	/*magneticField = new double[kgridNumber];
	for (int i = 0; i < kgridNumber;++i){
		magneticField[i] = 0;
	}*/
	detectedParticlesR1 = std::list<Particle*>();
	detectedParticlesR2 = std::list<Particle*>();

}

SpaceBin::~SpaceBin(){
	delete[] magneticField;
}

int SpaceBin::propagateParticle(Particle* particle ,double& time, double timeStep, const int rgridNumber){
	if(time != time){
		printf("aaa\n");
	}
	particle->setLocalMomentum(U);
	/*double lambda = getFreePath(particle);
	if(lambda != lambda){
		printf("aaa\n");
	}
	if(lambda == 0){
		printf("lambda == 0\n");
	}*/
	//double c2 = speed_of_light*speed_of_light;
	double gammaFactor = 1/sqrt(1 - U*U/c2);
	double colisionTime = getFreeTime(particle);
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

		double r = particle->absoluteX;

		return binByCoordinates(r,r1 - numberR*(r2 - r1), r2 - r1, rgridNumber);

	} else {
		largeAngleScattering(particle, time, timeStep);
		double particleR = particle->absoluteX;
		return binByCoordinates(particleR,0, r2 - r1, rgridNumber);
	}
}

void SpaceBin::largeAngleScattering(Particle* particle, double& time, double timeStep){
	double gammaFactor = 1/sqrt(1 - U*U/c2);
	//double lambda = getFreePath(particle);
	double colisionTime = getFreeTime(particle);
	double localV = particle->getLocalV();
	if(localV < epsilon){
		colisionTime = timeStep*2;
	}
	double r0;
	double sign;
	double absoluteV = particle->getAbsoluteVX();

	double deltat;
	
	if(abs(absoluteV) < epsilon){
		deltat = gammaFactor*(timeStep - time);
	}

	if(absoluteV >= 0){
		deltat = (r2 - particle->absoluteX)/absoluteV;
	} else {
		deltat = (r1 - particle->absoluteX)/absoluteV;
	}
	deltat = deltat*(1+epsilon);
	if(deltat < epsilon){
		deltat = epsilon;
	}
	if(deltat > gammaFactor*(timeStep - time)){
		deltat = gammaFactor*(timeStep - time);
	}

	time += deltat/gammaFactor;
	particle->absoluteX += particle->getAbsoluteVX()*deltat/gammaFactor;

	if(deltat < 0){
		printf("dt < 0\n");
		FILE* errFile = fopen("./output/errFile.dat","a");
		fprintf(errFile,"dt < 0\n");
		fprintf(errFile,"x= %lf, v= %lf, vx= %lf, r1= %lf, r2=%lf\n",particle->absoluteX, particle->getAbsoluteV(), absoluteV,r1,r2);
		fclose(errFile);
	}

	double probability;
	probability = 1 - exp(-deltat/colisionTime);
	if(uniRandom() <= probability){
		double cosLocalTheta = 2*(uniRandom() - 0.5);
		particle->localMomentumX = particle->localMomentum*cosLocalTheta;

		particle->setAbsoluteMomentum(U);
	}
}


int SpaceBin::binByCoordinates(double r, double r0, double deltar,const int rgridNumber){
	int i = lowerInt((r - r0)/deltar);

	if(i > rgridNumber - 1){
		//printf("aaa");
	}
	if( i < 0){
		//printf("aaa");
	}

	int result = i;
	return result;
}

/*double SpaceBin::getFreePath(Particle* particle){
	//return speed_of_light*particle.localMomentum/(particle.Z*electron_charge*B);
	double lambda = speed_of_light*particle->localMomentum/(particle->Z*electron_charge*B0);
	if( lambda != lambda){
		printf("aaa");
	}
	if( 0*lambda != 0*lambda){
		printf("aaa");
	}
	return lambda;
}*/

double SpaceBin::getFreeTime(Particle* particle){
	//return speed_of_light*particle.localMomentum/(particle.Z*electron_charge*B);
	double lambda = speed_of_light*particle->localMomentum/(particle->Z*electron_charge*B0);
	if( lambda != lambda){
		printf("aaa");
	}
	if( 0*lambda != 0*lambda){
		printf("aaa");
	}
	double time;
	//if(particle->localMomentum > particleLocalMomentum){
		//time = lambda*particle->mass/particleLocalMomentum;
	//} else {
	if(particle->getLocalV() < epsilon){
		time = defaultTimeStep*1000;
	} else {
		time = lambda*particle->mass/particle->getLocalV();
	}
	//}
	return time;
}

void SpaceBin::makeOneStep(Particle* particle, double colisionTime, double& time){
	double deltat = colisionTime/defaultTimeRelation;
	double particleU = particle->getAbsoluteV();
	double localMomentum = particle->localMomentum;
	double localMomentumX = particle->localMomentumX;

	double r = particle->absoluteX;


	double ux = particleU;

	//deltat = min4(deltat, abs(particle->absoluteX/ux), abs(particle->absoluteY/uy), abs(particle->absoluteZ/uz));
	if(deltat < epsilon){
		deltat = epsilon;
	}

	double tempx = particle->absoluteX + ux*deltat;


	double deltaR =  tempx - r;

	if(deltaR < 0){
		//printf("bbb");
	}


	while( (abs(deltaR) > (r2 - r1)/10) || (!(isInThisOrNear(tempx)) )){
		deltat = deltat/10;
		tempx = particle->absoluteX + ux*deltat;


		deltaR =  tempx - r;

		if(deltaR < 0){
			//printf("bbb");
		}
	}

	time += deltat;

	double maxTheta = sqrt(6*deltat/colisionTime);

	particle->absoluteX = tempx;


	scattering(particle, maxTheta);
}

void SpaceBin::scattering(Particle* particle, double maxTheta){
	double localTheta = acos(particle->localMomentumX/particle->localMomentum);
	if(abs(particle->localMomentum) < epsilon){
		localTheta = pi/2;
	}
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
	
	particle->localMomentumX = particle->localMomentum*cos(localTheta);
	particle->setAbsoluteMomentum(this);
}


void SpaceBin::resetDetectors(){
	detectedParticlesR1.clear();
	detectedParticlesR2.clear();


	particleMomentaZ.clear();
	particleWeights.clear();

	std::list<Particle*>::iterator it = particles.begin();
	while(it != particles.end()){
		Particle* particle = *it;
		delete particle;
		++it;
	}
	particles.clear();

	for(int j = 0; j < kgridNumber; ++j){
		//sortedParticles[j].clear();
	}

}

/////// I don't know how to do it






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

	return pmin;
}

void SpaceBin::sortParticles(double minK, double maxK){
	double kstep = (maxK - minK)/(kgridNumber - 1);
	std::list<Particle*>::iterator it = detectedParticlesR2.begin();
	while (it != detectedParticlesR2.end()){
		int i = (int) ( (**it).localMomentum - minK)/kstep;
		//sortedParticles[i].push_back(*it);
		++it;
	}
}

bool SpaceBin::isInBin(Particle* particle){
	double r = particle->absoluteX;
	if(r > r2) {
		return false;
	}
	if(r < r1) {
		return false;
	}
	return true;
}

bool SpaceBin::isInThisOrNear(double r){
	return true;
}

void SpaceBin::detectParticleR1(Particle* particle){

}

void SpaceBin::detectParticleR2(Particle* particle){

}

void SpaceBin::updateCosmicRayBoundMomentum(){

	if( particleMomentaZ.size() > 0){
		double maxp;
		maxp = particleMomentaZ[0];
		double weight = 0;
		weight += particleWeights[0];
		for(int i = 1; i < particleMomentaZ.size(); ++i){
			if(maxp < abs(particleMomentaZ[i])){
				maxp = abs(particleMomentaZ[i]);
			}
			weight += particleWeights[i];
		}
		double* distribution = new double[pgridNumber];
		for (int i = 0; i < pgridNumber; ++i){
			distribution[i] = 0;
		}
		int particleNumber = 0;
		double mass;
		double deltap = (2*maxp)/(pgridNumber);

		for(int i = 0; i < particleMomentaZ.size(); ++i){
			double p = particleMomentaZ[i];
			double particleWeight = particleWeights[i];
			if (abs(p) >= maxp){
				continue;
			}
			for(int j =0; j< pgridNumber; ++j){
				if (p <-maxp + deltap*(j + 1)){
					distribution[j] += particleWeight/(weight*deltap);
					break;
				}
			}
		}

		double a = 0;
		for(int i = 1; i < pgridNumber - 1; ++i){
			a = a + distribution[i]*deltap;
			distribution[i] = (distribution[i-1] + distribution[i] + distribution[i+1])/3;
		}

		updateTemperature(distribution, deltap);

		double maxDistribution = 0;
		for(int i = pgridNumber/2; i < pgridNumber-2; ++i){
			if(maxDistribution < (distribution[i])){
				maxDistribution = distribution[i];
			}
		}
		double scale = maxDistribution/(maxwell(0.0, massProton, temperature)) ;

		if(numberR == 10){
			FILE* outPDF = fopen("./output/zpdf1.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 100000000000000000000.0*(-maxp + i*deltap), distribution[i], scale*maxwell(-maxp + (i + 1/2)*deltap - centralMomentum, massProton, temperature));
			} 
			fclose(outPDF);
		}
		if(numberR == 20){
			FILE* outPDF = fopen("./output/zpdf2.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 100000000000000000000.0*(-maxp + i*deltap), distribution[i], scale*maxwell(-maxp + (i + 1/2)*deltap - centralMomentum, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 30){
			FILE* outPDF = fopen("./output/zpdf3.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 100000000000000000000.0*(-maxp + i*deltap - centralMomentum), scale*distribution[i], maxwell(-maxp + (i + 1/2)*deltap - centralMomentum, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 40){
			FILE* outPDF = fopen("./output/zpdf4.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 100000000000000000000.0*(-maxp + i*deltap), distribution[i], scale*maxwell(-maxp + (i + 1/2)*deltap - centralMomentum, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 50){
			FILE* outPDF = fopen("./output/zpdf5.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 100000000000000000000.0*(-maxp + i*deltap), distribution[i], scale*maxwell(-maxp + (i + 1/2)*deltap - centralMomentum, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 60){
			FILE* outPDF = fopen("./output/zpdf6.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 100000000000000000000.0*(-maxp + i*deltap), distribution[i], scale*maxwell(-maxp + (i + 1/2)*deltap - centralMomentum, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 70){
			FILE* outPDF = fopen("./output/zpdf7.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 100000000000000000000.0*(-maxp + i*deltap), distribution[i], scale*maxwell(-maxp + (i + 1/2)*deltap - centralMomentum, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 80){
			FILE* outPDF = fopen("./output/zpdf8.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 100000000000000000000.0*(-maxp + i*deltap), distribution[i], scale*maxwell(-maxp + (i + 1/2)*deltap - centralMomentum, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 90){
			FILE* outPDF = fopen("./output/zpdf9.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 100000000000000000000.0*(-maxp + i*deltap), distribution[i], scale*maxwell(-maxp + (i + 1/2)*deltap - centralMomentum, massProton, temperature));
			}
			fclose(outPDF);
		}
		if(numberR == 249){
			FILE* outPDF = fopen("./output/zpdf0.dat","w");
			for(int i = 0; i < pgridNumber; ++i){
				fprintf(outPDF,"%lf %lf %lf\n", 100000000000000000000.0*(-maxp + i*deltap), distribution[i], scale*maxwell(-maxp + (i + 1/2)*deltap - centralMomentum, massProton, temperature));
			}
			fclose(outPDF);
		}

		/*FILE* outPDF = fopen(fileName,"w");
		for(int i = 0; i < pgridNumber; ++i){
			fprintf(outPDF,"%lf %lf %lf\n", 100000000000000000000.0*(-maxp + i*deltap), distribution[i], scale*maxwell(-maxp + (i + 1/2)*deltap - centralMomentum, massProton, temperature));
		}
		fclose(outPDF);*/


		double j = 0;

		for(int i = 0; i < pgridNumber; ++i){
			distribution[i] -= scale*maxwell(-maxp + (i + 1/2)*deltap - centralMomentum, massProton, temperature);
			j += weight*((-maxp + (i + 1/2)*deltap)/massProton)*electron_charge*distribution[i]*deltap;
		}

		crFlux = j/volume;

		delete[] distribution;
	}
}

void SpaceBin::updateTemperature(double* distribution, double deltap){
	double sigma;
	int zeroP = pgridNumber/2;
	int sigmaIndex = zeroP;
	double maxDistribution = 0;
	int maxDistributionIndex = 1;
	for(int i = zeroP; i < pgridNumber-2; ++i){
		if(maxDistribution < (distribution[i])){
			maxDistribution = distribution[i];
			maxDistributionIndex = i;
		}
	}

	for(int i = maxDistributionIndex + 1;i < pgridNumber-2; ++i){
		if((distribution[i - 1] + distribution[i] + distribution[i + 1])/3< maxDistribution/exp(1.0)){
			sigmaIndex = i - 1;
			break;
		}
	}
	double d1 = distribution[sigmaIndex];
	double d2 = distribution[sigmaIndex + 1];


	//sigma =  (maxDistribution/exp(0.5) - d1)*deltap/(d2 - d1) + (sigmaIndex + 0.5 - zeroP)*deltap;
	sigma = (sigmaIndex + 0.5 - maxDistributionIndex)*deltap;
	temperature = sigma*sigma/(2*massProton*kBoltzman);
	temperature = temperature;
	centralMomentum = (maxDistributionIndex - zeroP)*deltap;

	//double minT = 0;
	//double maxT = (deltap*pgridNumber)*(deltap*pgridNumber)/(8*kBoltzman*massProton);
	//double T = findTemperature(minT,maxT,distribution, deltap);
	//temperature = T;
}


	


