#include "stdio.h"
#include "math.h"

#include "simulation.h"
#include "util.h"
#include "constants.h"
#include "specialmath.h"

void Simulation::filterHighHarmonics(){
	if(debugMode){
		printf("filtering harmonics\n");
	}
	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < znumber; ++j){
			for(int k = 0; k < znumber; ++k){
				double kx = 2*pi*i/(xsize);
				double ky = 2*pi*j/(ysize);
				double kz = 2*pi*k/(zsize);
				Vector3d kw = Vector3d(kx, ky, kz);
				if(kw.norm() > 0.3*2*pi/deltaX){
					if(debugMode){
						printf("harmonic %d %d %d\n", i, j, k);
					}
					Vector3d E = Vector3d(kx, -ky, 0);
					if(kx == 0 && ky == 0){
						E = Vector3d(0, 0, -kz);
					}
					E = E/E.norm();
					filterHarmonic(kw, E);
					E = kw.vectorMult(E);
					E = E/E.norm();
					filterHarmonic(kw, E);
				}
			}
		}
	}

	updateBoundaries();
}

void Simulation::filterHarmonic(Vector3d kw, Vector3d E){

	double imCoef = 0;
	double reCoef = 0;

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				imCoef += -EfieldX[i][j][k]*E.x*sin(kw.x*(xgrid[i] + deltaX/2) + kw.y*ygrid[j] + kw.z*zgrid[k]);
				reCoef += EfieldX[i][j][k]*E.x*cos(kw.x*(xgrid[i] + deltaX/2) + kw.y*ygrid[j] + kw.z*zgrid[k]);

				imCoef += -EfieldY[i][j][k]*E.y*sin(kw.x*xgrid[i] + kw.y*(ygrid[j] + deltaY/2) + kw.z*zgrid[k]);
				reCoef += EfieldY[i][j][k]*E.y*cos(kw.x*xgrid[i] + kw.y*(ygrid[j] + deltaY/2) + kw.z*zgrid[k]);

				imCoef += -EfieldZ[i][j][k]*E.z*sin(kw.x*xgrid[i] + kw.y*ygrid[j] + kw.z*(zgrid[k] + deltaZ/2));
				reCoef += EfieldZ[i][j][k]*E.z*cos(kw.x*xgrid[i] + kw.y*ygrid[j] + kw.z*(zgrid[k] + deltaZ/2));
			}
		}
	}

	double norm = cube(2*pi);

	imCoef /= norm;
	reCoef /= norm;

	for(int i = 0; i < xnumber; ++i){
		for(int j = 0; j < ynumber; ++j){
			for(int k = 0; k < znumber; ++k){
				EfieldX[i][j][k] += E.x*(imCoef*sin(kw.x*(xgrid[i] + deltaX/2) + kw.y*ygrid[j] + kw.z*zgrid[k]) - reCoef*cos(kw.x*(xgrid[i] + deltaX/2) + kw.y*ygrid[j] + kw.z*zgrid[k]));
				EfieldY[i][j][k] += E.y*(imCoef*sin(kw.x*xgrid[i] + kw.y*(ygrid[j] + deltaY/2) + kw.z*zgrid[k]) - reCoef*cos(kw.x*xgrid[i] + kw.y*(ygrid[j] + deltaY/2) + kw.z*zgrid[k]));
				EfieldZ[i][j][k] += E.z*(imCoef*sin(kw.x*xgrid[i] + kw.y*ygrid[j] + kw.z*(zgrid[k] + deltaZ/2)) - reCoef*cos(kw.x*xgrid[i] + kw.y*ygrid[j] + kw.z*(zgrid[k] + deltaZ/2)));
			}
		}
	}
}