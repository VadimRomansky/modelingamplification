#include "stdio.h"
#include <stdlib.h>
#include <time.h>   

#include "constants.h"
#include "matrix3d.h"
#include "util.h"
#include "simulation.h"

//test solving matrix

/*int main(){
	double** leftHalf;
	double* rightPart;
	double* output;

	leftHalf = new double*[6];
	rightPart = new double[6];
	output = new double[6];

	for(int i = 0; i < 6; ++i){
		leftHalf[i] = new double[3];
		for(int j = 0; j < 3; ++j){
			leftHalf[i][j] = 0;
		}
		rightPart[i] = 1;
		output[i] = 0;
	}
	
	
	//	1 2 3 0 0 0      1
	//	2 1 3 0 0 0      1
	//	4 1 1 0 0 0      1
	//  0 3 0 1 0 0      1
	//	0 0 2 0 1 0      1
	//	2 0 0 0 0 0      1

	//	x = 
	
	//	0.16667
	//	0.16667
	//	0.16667
	//	0.5
	//	0.66667
	//	0.66667

	leftHalf[0][0] = 1;
	leftHalf[1][1] = 1;
	leftHalf[2][2] = 1;
	leftHalf[5][0] = 2;
	leftHalf[3][1] = 3;
	leftHalf[4][2] = 2;
	leftHalf[0][1] = 2;
	leftHalf[0][2] = 3;
	leftHalf[1][0] = 2;
	leftHalf[1][2] = 3;
	leftHalf[2][0] = 4;
	leftHalf[2][1] = 1;


	solveSpecialMatrix(leftHalf, rightPart, output);

	for(int i = 0; i < 6; ++i){
		printf("%lf\n", output[i]);
	}

	for(int i = 0; i < 6; ++i){
		delete[] leftHalf[i];
	}
	delete[] leftHalf;
	delete[] rightPart;
	delete[] output;

	return 0;
}*/


//test McDonald
int main()
{
	srand (time(NULL));
	//double x = uniformDistribution()*10;
	double x = 1;
	double index = 2;

	double result = McDonaldFunction(x, index);

	printf("K(%lf, %lf) = %g\n", index, x, result);
}

/*int main()
{	
	srand (time(NULL));

	Simulation simulation;

	simulation.simulate();

}*/

