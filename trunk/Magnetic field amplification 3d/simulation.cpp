#include "stdafx.h"
#include "simulation.h"
#include "SpaceBin.h"
#include "particle.h"
#include "output.h"
#include "util.h"

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
	/*for(int i = 0; i < xgridNumber; ++i){
		delete[] pressureSpectralDensity[i];
	}
	delete[] pressureSpectralDensity;*/
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
					u = U0*sqr((upstreamR + deltaR/2)/R);
					//density = density0*sqr((upstreamR + deltaR/2)/R);
				} else {
					u = U0*sqr((upstreamR + deltaR/2)/R)/Rtot;
					//density = density0*sqr((upstreamR + deltaR/2)/R)/Rtot;
				}
				bins[i][j][k] = new SpaceBin(R,Theta,Phi,deltaR,deltaTheta,deltaPhi,u,density,Theta,Phi,temperature,B0,i,j,k);
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
///// главная функция программы.
void Simulation::simulate(){
	FILE* outIteration = fopen("./output/tamc_iteration.dat","w");
	fclose(outIteration);
	srand ( time(NULL) );
	initializeProfile();
	//FILE* outIteration = fopen("./output/tamc_iteration.dat","w");
	introducedParticles = getParticles();
	for (int itNumber  = 0; itNumber < iterationNumber; ++itNumber){ 
		printf("%s", "\n");
		int j = 0;
		int l = 0;
		printf("%s", "Iteration started\n");
		printf("%s", "Particle propagation\n");
		std::list<Particle*>::iterator it = introducedParticles.begin();
		while( it != introducedParticles.end()){
			Particle* particle = *it;
			printf("%d %s",l,"\n");
			++l;
			bool side = false;
			double r = sqrt(particle->absoluteX*particle->absoluteX + particle->absoluteY*particle->absoluteY + particle->absoluteZ*particle->absoluteZ);
			double theta = acos(particle->absoluteZ/r);
			double phi =atan2(particle->absoluteY, particle->absoluteX);
			if(phi < 0){
				phi = phi + 2*pi;
			}
			int* index = SpaceBin::binByCoordinates(r, theta, phi,upstreamR,deltaR,deltaTheta,deltaPhi);
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
				int* tempIndex = bin->propagateParticle(particle,time, timeStep);
				delete[] index;
				index = tempIndex;
				if (time >= timeStep){
					break;
				}
			}
			++it;
		}
		/*printf("%s", "Fluxes updating\n");
		for (int i = 0; i < rgridNumber; ++i){
			for (int j = 0; j < phigridNumber; ++j){
				for (int k = 0; k < thetagridNumber; ++k){
					bins[i][j][k]->updateFluxes();
				}
			}
		}*/
		printf("%s", "Reseting profile\n");
		resetProfile();
		printf("%s", "magnetic Field updating\n");
		//updateMagneticField();
		output(*this);
		resetDetectors();
		printf("%s","iteration № ");
		printf("%d\n",itNumber);
		outIteration = fopen("./output/tamc_iteration.dat","a");
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
		outputPDF(list,"./output/tamc_pdf_down.dat",*this,minp,maxp);
		outputStartPDF(list,"./output/tamc_pdf_start.dat",*this,minp,maxp);
		outputRadialProfile(bins,0,0);
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
		for(int i = 1; i < rgridNumber; i++){
			double massFlux1 = bins[i - 1][0][0]->particleMassFlux.fluxR2;
			double massFlux2 = 0;
			if ( i < rgridNumber){
				massFlux2 = bins[i + 1][0][0]->particleMassFlux.fluxR1;
			}
			double deltaM = massFlux2 + massFlux1 - bins[i][0][0]->particleMassFlux.fluxR1 - bins[i][0][0]->particleMassFlux.fluxR2;
			bins[i][0][0]->density += deltaM/(bins[i][0][0]->volume);
		}
	}
}

////обнуляет счётчики зарегистрированных частиц
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

/*std::list <Particle> Simulation::getParticleUniformDistribution(double  pmax,int  pnumber, int thetanumber){
	double pstep = pmax/pnumber;
	double thetastep = pi/(thetanumber);
	std::list <Particle> l = std::list <Particle>();
	for(int i=0; i<pnumber; ++i){
		for(int j=0; j<thetanumber; ++j){
			Particle particle = Particle(downstreamR, temperature,A,Z,U0,i*pstep+pstep/2,j*thetastep+thetastep/2);
			particle.weight = maxwell(i*pstep+pstep/2, temperature, A*massProton)*thetastep*pstep*sin(j*thetastep+thetastep/2)/(4*pi);
			l.push_front(particle);
		}
	}
	return l;
}*/
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
	/*updateMaxMinK();
	updatePressureSpectralDensity();
	for(int i =0; i < kgridNumber; ++i){
		double k = minK + (maxK-minK)*i/(kgridNumber-1);
		for(int j = 0; j < xgridNumber; ++j){
			xbins[j]->magneticField[i] = turbulenceSeed*turbulenceSeed/(4*pi*maxK*k*log(maxK/minK));
		}
	}
	Xbin* bin1 = NULL;
	Xbin* bin2 = xbins[0];
	bin2->updateMagneticField();
	for(int i = 1; i <xgridNumber; ++i){
		bin1 = bin2;
		bin2 = xbins[i];
	    double dx = (bin2->right - bin1->left)/2;
		double gradU = gradientSpeed(i);
		//printf("%s %d","bin number",i);
		//printf("\n");
		evaluateMagneticField(bin1->magneticField,bin2->magneticField,dx,1000,gradU,bin1->density,bin1->U,i);
		bin2->updateMagneticField();
	}*/
}

////градиент скорости
vector3d Simulation::gradientSpeed(int i, int j, int k){
	/*Xbin* bin1;
	Xbin* bin2;
	if(i == 0){
		bin1 = xbins[1];
		bin2 = xbins[0];
	} else {
		bin1 = xbins[i];
		bin2 = xbins[i-1];
	}
	///// расстояние между центрами ящиков
	return (bin1->U-bin2->U)*2/(bin1->right-bin2->left);*/
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
	/*Xbin* bin;
	double instability = 0;
	bin = xbins[i];
	double k = minK + j*(maxK - minK)/(kgridNumber-1);
	double va = B0/sqrt(4*pi*bin->density);
	if (resonantInstability){
		/////TODO Z?		
		double dx = (xbins[i]->right - xbins[i-1]->left)/2;
		double gradp = 0;
		//double p1 = xbins[i]->pressureSpectralDensity(k, B0, density0, U0, Z, A);
		//double p2 = xbins[i-1]->pressureSpectralDensity(k, B0, density0, U0, Z, A);
		double p1 = pressureSpectralDensity[i][j];
		double p2 = pressureSpectralDensity[i-1][j];
		gradp = (p1 -p2)/dx;
		/////TODO B0?
		double pres = Z*electron_charge*B0/(speed_of_light*k);
		/////TODO как считать давление космических лучей?
		instability += va*gradp*pres/(k);
		if(instability != instability){
			printf("NaN instability");
		}
	}
	if (bellInstability){
		/////TODO как обработать нижнюю границу по k?
		double j = bin->crMassFlux*electron_charge*Z/(A*massProton);
		double kc = B0*j/(speed_of_light*bin->density*va*va);
		if(k < kc){
			instability += 2*va*k*sqrt(kc/k - 1.0)*w;
			if(instability != instability){
				printf("NaN instability");
			}
		}
	}
	double x = instability*0;
	if((instability != instability)||(x != x)){
		printf("NaN instability");
	}
	return instability;*/
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
	/*double va;
	double kd;
	Xbin* bin = xbins[i];
	//////TODO Z и A?
	kd = Z*electron_charge*B0/(speed_of_light*sqrt(massProton*A*kBoltzman*bin->temperature));
	//////TODO B0 ?
	va = B0/sqrt(4*pi*bin->density);
	return va*k*k*w/kd;*/
	return 0;
}


void Simulation::updateMaxMinK(){
	/*double pmin;
	double pmax;
	double p1;
	double p2;
	pmax = xbins[0]->getMaxLocalP();
	pmin = xbins[0]->getMinLocalP();
	for(int i = 1; i < xgridNumber; ++i){
		p1 = xbins[i]->getMaxLocalP();
		p2 = xbins[i]->getMinLocalP();
		if(p1 > pmax){
			pmax = p1;
		}
		if(p2 < pmin){
			pmin = p2;
		}
	}
	minK = Z*electron_charge*B0/(speed_of_light*pmax);
	maxK = Z*electron_charge*B0/(speed_of_light*pmin);
	minK = minK/1.1;
	maxK = maxK*1.1;
	for(int i = 0; i < xgridNumber; ++i){
		xbins[i]->minK = minK;
		xbins[i]->maxK = maxK;
	}*/
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

void Simulation::evaluateMagneticField(double* startField,double* endField,double deltax,int gridNumber,double gradU, double density,double U,int binNumber){
	/*double* currentField= new double[kgridNumber];
	double* newField= new double[kgridNumber];
	for(int i = 0; i<kgridNumber; ++i){
		currentField[i] = startField[i];
	}
	double dx = deltax/gridNumber;
	for(int i = 0; i<gridNumber; ++i){
		if (i%100 ==0){
			//printf("%d %s",i,"\n");
		}
		for(int j = 0; j < kgridNumber; ++j){
			double w = currentField[j];
			double k = minK + j*(maxK-minK)/(kgridNumber-1);
			double rho = density;
			////TODO i-1?
			double x = w*0;
			if ((w != w)||(x != x)){
				printf("%s","NaN w");
			}
			double cDW = cascadingDerivativeW(w,k,rho);
			x = cDW*0;
			if ((cDW != cDW)||(x != x)){
				printf("%s", "NaN cDW\n");
			}
			double dFK = derivativeFieldK(currentField,j);
			x = dFK*0;
			if ((dFK != dFK)||(x != x)){
				printf("%s", "NaN dFK\n");
			}
			double inst = instability(w,binNumber,j);
			x = inst*0;
			if ((inst != inst)||(x != x)){
				printf("%s", "NaN inst\n");
			}
			double cDK =  cascadingDerivativeK(w,k,rho);
			x = cDK*0;
			if ((cDK != cDK)||(x != x)){
				printf("%s", "NaN cDK\n");
			}
			double diss = dissipation(w,binNumber,k);
			x = diss*0;
			if ((diss != diss)||(x != x)){
				printf("%s", "NaN diss\n");
			}
			if( w != w){
				printf("%s", "NaN w\n");
			}
			if( gradU != gradU){
				printf("%s","NaN gradU\n");
			}
			//double dw = (dx/bin1->U)*((delta*cascadingDerivativeW(w,k,rho)-beta*k*gradU)*derivativeFieldK(i-1,j)-(alpha-beta)*w*gradU + gamma*instability(w,i,j) -delta*cascadingDerivativeK(w,k,rho) - epsilon*dissipation(w,i,k));
			double dw = (dx/U)*((-delta*cDW+beta*k*gradU)*dFK-(alpha-beta)*w*gradU + gamma*inst -delta*cDK - epsilon*diss);
			x = w*0;
			if(( w != w)||(x != x)){
				printf("%s", "NaN w\n");
			}
			x = dw*0;
			if(( dw != dw)||(x != x)){
				printf("%s", "NaN w\n");
			}
			if((w < 1E100)&&(dw < 1E100)){
				newField[j] = w + dw;
			} else {
				newField[j] = w;
				//printf("%s","w goes to infinity");
			}
			x = newField[j]*0;
			if((newField[j] != newField[j])||(x != x)){
				printf("%s","NaN newField");
			}
			if(newField[j] < 0){
				newField[j] = 0;
				//printf("%s","field < 0");
			}
	
		}

		for(int j = 0; j<kgridNumber; ++j){
			currentField[j] = newField[j];
		}
	}
	for(int i = 0; i<kgridNumber; ++i){
		endField[i] = currentField[i];
	}*/
}

void Simulation::updatePressureSpectralDensity(){
	/*for(int i = 0; i < xgridNumber; ++i){
		xbins[i]->sortParticles(minK,maxK);
		//printf("%d %s",i,"\n");
		for(int j = 0; j < kgridNumber; ++j){
			double k = minK + j*(maxK - minK)/(kgridNumber-1);
			pressureSpectralDensity[i][j] = xbins[i]->pressureSpectralDensity(j, density0, U0, Z, A); 
		}
	}*/
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
				for( int l = 0; l < particlesNumber; ++l){
					//if(order(bin->phi1, phi, bin->phi2) && order(bin->theta1, theta, bin->theta2) && order(bin->r1, r, bin->r2)){
						++allParticlesNumber;
						printf("%d",allParticlesNumber);
						printf("%s","\n");
						Particle* particle = new Particle( A, Z,bin, true);
						particle->weight /= particlesNumber;
						startPDF.push_front(particle);
						list.push_front(particle);
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
	int j = trunc(theta / deltaTheta);
	int k = trunc(phi / deltaPhi);
	return bins[0][j][k];
}

void Simulation::detectFromTo(int fromIndexR, int fromIndexTheta, int fromIndexPhi, int toIndexR, int toIndexTheta, int toIndexPhi, const Particle& particle){
	/*if(fromIndexR > toIndexR){
		if (fromIndexR < rgridNumber) bins[fromIndexR][fromIndexTheta][fromIndexPhi]->detectedParticlesR1.push_back(Particle(particle));
		if (toIndexR >=0) bins[toIndexR][toIndexTheta][toIndexPhi]->detectedParticlesR2.push_back(Particle(particle));
		return;
	}
	if(fromIndexR < toIndexR){
		if (fromIndexR >= 0) bins[fromIndexR][fromIndexTheta][fromIndexPhi]->detectedParticlesR2.push_back(Particle(particle));
		if (toIndexR < rgridNumber) bins[toIndexR][toIndexTheta][toIndexPhi]->detectedParticlesR1.push_back(Particle(particle));
		return;
	}
	if(fromIndexTheta != toIndexTheta){
		if(fromIndexTheta == 0){
			if(toIndexTheta == thetagridNumber - 1){
			} else {
				bins[fromIndexR][fromIndexTheta][fromIndexPhi]->detectedParticlesTheta2.push_back(Particle(particle));	
			}
		}
	}
	if(fromIndexPhi != toIndexPhi){
	}
	if(fromIndexTheta 
	if*/
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
