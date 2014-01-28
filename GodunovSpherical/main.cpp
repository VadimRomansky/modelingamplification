#include "stdafx.h"
//#include <crtdbg.h>
#include "input.h"
#include "simulation.h"
#include "input.h"
#include <omp.h>
#include <time.h>


int main()
{
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_CHECK_ALWAYS_DF);
	//omp_set_num_threads( 12 );
	Simulation* simulation = readInput();
	simulation->simulate();
	return 0;
}
