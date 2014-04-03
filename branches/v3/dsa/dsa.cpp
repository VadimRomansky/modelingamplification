#include "math.h"
#include "stdio.h"
#include "solver.h"

int main(){

	int Nx=200, Np=100, Nt = 300000;
	double   a=10000, b=10000, Pmin=0.01, Pmax=100000;
	double** gn = new double*[Nx];
	double** g = new double*[Nx];
	for(int i = 0; i < Nx; ++i){
		gn[i] = new double[Np-1];
		g[i] = new double[Np-1];
	}
	double* x = new double[Nx];
	double dt, h1, h2, ymin, dy, y;


	ymin=log(Pmin);
	dy=log(Pmax/Pmin)/Np;

	//создание сетки по х
	h1=0.5*Nx/log(1.0+a);
	h2=0.5*Nx/log(1.0+b);
	for(int i=0; i < Nx/2; ++ i){ 
		x[i] = 1 - exp(-(1.0*(i+1)-0.5*Nx)/h1);
	}
	for(int i=Nx/2; i < Nx; ++i){
		x[i] = exp((1.0*(i+1)-0.5*Nx)/h1)-1.0;
	}

	dt=0.001;
	printf("%lf %g %g\n",dt, x[Nx/2-1]-x[Nx/2-2], x[Nx/2]-x[Nx/2-1]);

	for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Np-1; ++j){
			gn[i][j] = 0;
			g[i][j] = 0;
		}
	}
	int iter = 0;
	
	for(iter = 1; iter <= Nt; ++iter){
		iter = iter + 1;
		solver(a, ymin, x, dy, dt, Nx, Np, gn, g);
		for(int i = 0; i < Nx; ++i){
			for(int j = 0; j <Np; ++j){
				gn[i][j] = g[i][j];
			}
		}
		printf("%d\n",iter);
	}


	FILE* fp4 = fopen("data/fp4.dat","w");
	for(int k = 0; k < Np-1; ++k){
		y=ymin+k*dy;
		fprintf(fp4,"%d %g %g\n",k, exp(y), g[Nx/2-1][k]*exp(y));
	}
	fclose(fp4);

	FILE* NXFile = fopen("data/Nx.dat","w");
	fprintf(NXFile, "%d", Nx);
	fclose(NXFile);

	FILE* NPFile = fopen("data/Np.dat","w");
	fprintf(NPFile, "%d", Np);
	fclose(NPFile);

	FILE* gfull = fopen("data/gfull.dat","w");
	for(int i = 0; i < Nx; ++i){
		for(int k = 0; k < Np-1; ++k){
			fprintf(gfull, "%lf\n", g[i][k]);
		}
	}
	fclose(gfull);

	FILE* xfile = fopen("data/x.dat","w");
	for(int i = 0; i < Nx; ++i){
		fprintf(xfile, "%lf\n", x[i]);
	}
	fclose(xfile);
}