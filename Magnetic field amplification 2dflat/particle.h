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

	double absoluteMomentum;
	double absoluteMomentumX;
	double absoluteMomentumY;

	double localMomentum;
	double localMomentumX;
	double localMomentumY;

	int number;


	bool  isCosmicRay;
	double initialMomentum;
	double initialLocalMomentum;
	double previousAbsoluteMomentum;
	double weight;
	std::list <double> path;
	bool writePath;
	Particle();
	Particle(const Particle& p);
	Particle(int a, int znumber, SpaceBin* bin, bool wpath, int n); 
	double getAbsoluteV();
	double getAbsoluteVX();
	double getAbsoluteVY();
	void setAbsoluteMomentum(double Ux, double Uy);
	void setLocalMomentum(double Ux, double Uy);
	void setAbsoluteMomentum(SpaceBin* bin);
	void setLocalMomentum(SpaceBin* bin);
	double getLocalV();
	double getLocalVX();
	double getLocalVY();
	double getEnergy();
	double getInitialEnergy();

	void moveToBinLeftX(SpaceBin* bin);
	void moveToBinLeftY(SpaceBin* bin);
};
#endif