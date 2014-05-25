// TestOMP.cpp : Defines the entry point for the console application.
//

#include "omp.h"
#include "stdio.h"


int main()
{
	omp_set_num_threads(2);

	double k=1;
	double* a = new double[1000];
	for(int i = 0; i < 1000; ++i){
		a[i] = 0;
	}
	double t1 = omp_get_wtime();
	#pragma omp parallel for schedule(static,500)
	for(int i = 0; i < 1000; ++i){
		for(int j = 0; j < 10000; ++j){
			a[i] += j;
		}
	}
	double t2 = omp_get_wtime();
	printf("parllel %lf\n",t2-t1);
	printf("%lf\n", k);
	k=1;
	for(int i = 0; i < 1000; ++i){
		a[i] = 0;
	}
	t1 = omp_get_wtime();
	for(int i = 0; i < 1000; ++i){
		for(int j = 0; j < 10000; ++j){
			a[i] += j;
		}
	}
	t2 = omp_get_wtime();
	printf("not parllel %lf\n",t2-t1);
	printf("%lf\n", k);

	return 0;
}

