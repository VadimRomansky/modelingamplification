#ifndef PARTICLE_H
#define PARTICLE_H

#include <list>
#include "stdafx.h"
#include "constants.h"

class SpaceBin;

class Particle{
public:
	int A;
	int Z;
	double mass;
	double absoluteX;
	double absoluteY;
	double absoluteZ;

	double absoluteMomentum;
	double absoluteMomentumTheta;
	double absoluteMomentumPhi;

	double localMomentum;
	double localMomentumX;
	double localMomentumY;
	double localMomentumZ;


	bool  isCosmicRay;
	double initialMomentum;
	double initialLocalMomentum;
	double weight;
	std::list <double> path;
	bool writePath;
	Particle();
	Particle(const Particle& p);
	Particle(double r, double temperature, int a, int z);
	Particle(double x, double y, double z, double temperature, int a, int znumber);
	Particle(double x, double y, double z, double temperature, int a, int znumber, double Ur, double Utheta, double Uphi);
	Particle(double x, double y, double z, double temperature, int a, int znumber, double Ur, double Utheta, double Uphi, bool wPath);
	Particle(double x, double y, double z, double temperature, int a, int znumber, double Ur, double Utheta, double Uphi, double px, double py, double pz, bool wPath);
	Particle(double x, double y, double z, double temperature, int a, int znumber, double Ur, double Utheta, double Uphi, double px,double py, double pz, double ptheta, double pphi);
	Particle(int a, int znumber, SpaceBin* bin, bool wpath); 
	void LorentzTransition(double v1, double v2);
	double getAbsoluteV();
	void setAbsoluteMomentum(double U, double Utheta, double Uphi);
	void setLocalMomentum(double U, double Utheta,double Uphi);
	void setAbsoluteMomentum(SpaceBin* bin);
	void setLocalMomentum(SpaceBin* bin);
	double getLocalV();
	double getAbsoluteR();
	double getAbsoluteTheta();
	double getAbsolutePhi();
	double getEnergy();
	double getRadialSpeed();
	void moveToBinRight(SpaceBin* bin);
};
#endif