#include "constants.h"
#include "matrix3d.h"
#include "stdio.h"

Matrix3d sum(Matrix3d& m1, Matrix3d& m2){
	Matrix3d newMatrix;
	for(int i = 0; i < 3; ++i){
		for(int j = 0; j < 3; ++j){
			newMatrix.matrix[i][j] = m1.matrix[i][j] + m2.matrix[i][j];
		}
	}
	return newMatrix;
}

int main()
{	
	Matrix3d A;
	Matrix3d B;
	Matrix3d C;
	for(int i = 0; i < 10000000; ++i){
		A = B+C;
		if(i%1000 == 0){
			printf("%d\n",i);
		}
	}
	return 0;
}

