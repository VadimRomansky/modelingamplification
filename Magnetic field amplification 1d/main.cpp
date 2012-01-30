// Magnetic field amplification.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
//#include <crtdbg.h>
#include "input.h"
#include "simulation.h"




int main()
{
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_CHECK_ALWAYS_DF);
	srand ( time(NULL) );
	Simulation* simulation = readInput();
	simulation->simulate();
	return 0;
}

