#include "stdio.h"
#include "simulation.h"

#include "input.h"

Simulation readInput(FILE* inputFile) {
	int xnumber;
	char ch = ' ';
	fscanf(inputFile, "%d", &xnumber);
	int ynumber;
	fscanf(inputFile, "%d", &ynumber);
	int znumber;
	fscanf(inputFile, "%d", &znumber);

	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double xsize;
	fscanf(inputFile, "%lf", &xsize);
	double ysize;
	fscanf(inputFile, "%lf", &ysize);
	double zsize;
	fscanf(inputFile, "%lf", &zsize);

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double temperature;
	fscanf(inputFile, "%lf", &temperature);

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double density;
	fscanf(inputFile, "%lf", &density);

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double Ex;
	fscanf(inputFile, "%lf", &Ex);

	double Ey;
	fscanf(inputFile, "%lf", &Ey);

	double Ez;
	fscanf(inputFile, "%lf", &Ez);

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double Bx;
	fscanf(inputFile, "%lf", &Bx);

	double By;
	fscanf(inputFile, "%lf", &By);

	double Bz;
	fscanf(inputFile, "%lf", &Bz);

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	int maxIterations;
	fscanf(inputFile, "%d", &maxIterations);

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	double maxTime;
	fscanf(inputFile, "%lf", &maxTime);

	ch = ' ';
	while(ch != '\n') {
		fscanf(inputFile, "%c", &ch);
	}

	int particlesPerBin;
	fscanf(inputFile, "%d", &particlesPerBin);

	return Simulation(xnumber, ynumber, znumber, xsize, ysize, zsize, temperature, density, Ex, Ey, Ez, Bx, By, Bz, maxIterations, maxTime, particlesPerBin);
}