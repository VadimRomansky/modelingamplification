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

	double absoluteMomentum;
	double absoluteMomentumX;

	double localMomentum;
	double localMomentumX;

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
	Particle(double x, double temperature, int a, int znumber);
	Particle(double x, double temperature, int a, int znumber, double Ur);
	Particle(double x, double temperature, int a, int znumber, double Ur, bool wPath);
	Particle(double x, double temperature, int a, int znumber, double Ur, double px, double py, double pz, bool wPath);
	Particle(double x, double temperature, int a, int znumber, double Ur, double px,double py, double p);
	Particle(int a, int znumber, SpaceBin* bin, bool wpath, int n); 
	void LorentzTransition(double v1, double v2);
	double getAbsoluteV();
	double getAbsoluteVX();
	void setAbsoluteMomentum(double U);
	void setLocalMomentum(double U);
	void setAbsoluteMomentum(SpaceBin* bin);
	void setLocalMomentum(SpaceBin* bin);
	double getLocalV();
	double getLocalVX();
	double getEnergy();
	double getInitialEnergy();
	void moveToBinRight(SpaceBin* bin);
	void moveToBinLeft(SpaceBin* bin);
};
#endif