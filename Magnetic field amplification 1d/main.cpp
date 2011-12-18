// Magnetic field amplification.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "input.h"
#include "simulation.h"




int main()
{
	srand ( time(NULL) );
	Simulation* simulation = readInput();
	simulation->simulate();
	return 0;
}

