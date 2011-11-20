// Magnetic field amplification.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "input.h"
#include "simulation.h"




int main()
{
	//printf("Hello, world!");
	/*const Constant::massProton = 1;*/
	srand ( time(NULL) );
	Simulation* simulation = readInput();
	simulation->simulate();
	/*FILE* file = fopen("test.txt","w");
	double* a = new double[100];
	for(int i = 0 ; i <100; ++i){
		a[i] = 0;
	}
    for(int i = 0; i < 100000; ++i){
		double d = uniRandom();
		int j = d*100;
		a[j] = a[j] +1.0/10000;
	}
	for(int i = 0; i < 100; ++i){
		fprintf(file,"%lf %s",a[i],"\n");
	}
	fclose(file);*/
	return 0;
}

