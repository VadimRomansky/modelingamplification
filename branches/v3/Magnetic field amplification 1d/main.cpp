// Magnetic field amplification.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
//#include <crtdbg.h>
#include "input.h"
#include "simulation.h"
#include <omp.h>




int main()
{
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_CHECK_ALWAYS_DF);
	omp_set_num_threads( 2 );
	srand ( time(NULL) );
	Simulation* simulation = readInput();
	simulation->simulate();
	return 0;
}

