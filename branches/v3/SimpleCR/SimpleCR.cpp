// SimpleCR.cpp : Defines the entry point for the console application.
//
#include "simulation.h"
#include "input.h"



int main()
{
	Simulation* simulation = readInput();
	simulation->simulate();
	return 0;
}

