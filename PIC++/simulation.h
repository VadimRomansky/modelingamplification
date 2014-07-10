#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "stdlib.h"
#include "list"

class Particle;

class Simulation{
	int xnumber;
	int ynumber;
	int znumber;

	double*** Efield;
	double*** Bfield;

	std::list<Particle*> particles;

	Simulation();
};

#endif