#include "stdlib.h"
#include "stdio.h"

#include "simulation.h"
#include "util.h"
#include "particle.h"
#include "constants.h"
#include "matrix3d.h"
#include "specialmath.h"

void Simulation::evaluateFields() {
	printf("evaluating fields\n");

	updateBfield(deltaT/2);

	updateEfield(deltaT);

	updateBfield(deltaT/2);

	updateBoundaries();

	double error = 4*pi*chargeDensity[xnumber/2][ynumber/2][znumber/2] - evaluateDivE(xnumber/2, ynumber/2, znumber/2);
	cleanupDivergence();
	error = 4*pi*chargeDensity[xnumber/2][ynumber/2][znumber/2] - evaluateDivE(xnumber/2, ynumber/2, znumber/2);
}

void Simulation::updateEfield(double dt) {
	for(int i = 0; i < xnumber; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				updateEfieldX(i, j, k, dt);
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				updateEfieldY(i, j, k, dt);
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber; ++k) {
				updateEfieldZ(i, j, k, dt);
			}
		}
	}
}

void Simulation::updateEfieldX(int i, int j, int k, double dt) {
	int prevJ = j - 1;
	if(prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int prevK = k - 1;
	if(prevK < 0) {
		prevK = znumber - 1;
	}
	if(j >= ynumber) {
		j = 0;
	}
	if(k >= znumber) {
		k = 0;
	}

	int middleJ = j;
	if(middleJ >= ynumber) {
		middleJ = 0;
	}
	int middleK = k;
	if(middleK >= znumber) {
		middleK = 0;
	}

	double BrightY = BfieldZ[i][middleJ][middleK];
	double BleftY = BfieldZ[i][prevJ][middleK];

	double BrightZ = BfieldY[i][middleJ][middleK];
	double BleftZ = BfieldY[i][middleJ][prevK];

	double rotBx = (BrightY - BleftY)/deltaY - (BrightZ - BleftZ)/deltaZ;

	EfieldX[i][j][k] += (speed_of_light_normalized*rotBx - 4*pi*electricFluxX[i][j][k])*dt;

	alertNaNOrInfinity(EfieldX[i][j][k], "EfieldX = NaN\n");
}

void Simulation::updateEfieldY(int i, int j, int k, double dt) {
	int middleI = i;
	int prevI = i-1;
	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		if(i == 0) {
			EfieldY[i][j][k] = 0;
			return;
		}
		if(i == xnumber) {
			EfieldY[i][j][k] = E0.y;
			return;
		}
	} else if(boundaryConditionType == PERIODIC) {
		if(middleI >= xnumber) {
			middleI = 0;
		}
		if(prevI < 0) {
			prevI = xnumber - 1;
		}
	}

	int prevJ = j - 1;
	if(prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int prevK = k - 1;
	if(prevK < 0) {
		prevK = znumber - 1;
	}
	int middleJ = j;
	if(middleJ >= ynumber) {
		middleJ = 0;
	}
	int middleK = k;
	if(middleK >= znumber) {
		middleK = 0;
	}

	double BrightX = BfieldZ[middleI][middleJ][middleK];
	double BleftX = BfieldZ[prevI][middleJ][middleK];

	double BrightZ = BfieldX[middleI][middleJ][middleK];
	double BleftZ = BfieldX[middleI][middleJ][prevK];

	double rotBy = -(BrightX - BleftX)/deltaX + (BrightZ - BleftZ)/deltaZ;

	EfieldY[i][j][k] += (speed_of_light_normalized*rotBy - 4*pi*electricFluxY[i][j][k])*dt;

	alertNaNOrInfinity(EfieldY[i][j][k], "EfieldY = NaN\n");
}

void Simulation::updateEfieldZ(int i, int j, int k, double dt) {
	int middleI = i;
	int prevI = i-1;
	if(boundaryConditionType== SUPERCONDUCTERLEFT){
		if(i == 0) {
			EfieldZ[i][j][k] = 0;
			return;
		}
		if(i == xnumber) {
			EfieldZ[i][j][k] = E0.z;
			return;
		}
	} else if(boundaryConditionType == PERIODIC) {
		if(middleI >= xnumber) {
			middleI = 0;
		}
		if(prevI < 0) {
			prevI = xnumber - 1;
		}
	}

	int prevJ = j - 1;
	if(prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int prevK = k - 1;
	if(prevK < 0) {
		prevK = znumber - 1;
	}
	if(j >= ynumber) {
		j = 0;
	}
	if(k >= znumber) {
		k = 0;
	}

	int middleJ = j;
	if(middleJ >= ynumber) {
		middleJ = 0;
	}
	int middleK = k;
	if(middleK >= znumber) {
		middleK = 0;
	}

	double BrightX = BfieldY[middleI][middleJ][middleK];
	double BleftX = BfieldY[prevI][middleJ][middleK];

	double BrightY = BfieldX[middleI][middleJ][middleK];
	double BleftY = BfieldX[middleI][prevJ][middleK];

	double rotBy = (BrightX - BleftX)/deltaX - (BrightY - BleftY)/deltaY;

	EfieldZ[i][j][k] += (speed_of_light_normalized*rotBy - 4*pi*electricFluxY[i][j][k])*dt;

	alertNaNOrInfinity(EfieldZ[i][j][k], "EfieldZ = NaN\n");
}

void Simulation::updateBfield(double dt) {
	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				updateBfieldX(i, j, k, dt);
			}
		}
	}

	for(int i = 0; i< xnumber; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			for(int k = 0; k < znumber; ++k) {
				updateBfieldY(i, j, k, dt);
			}
		}
	}

	for(int i = 0; i< xnumber; ++i) {
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber + 1; ++k) {
				updateBfieldZ(i, j, k, dt);
			}
		}
	}	
}


void Simulation::updateBfieldX(int i, int j, int k, double dt) {
	int nextJ = j + 1;
	if(nextJ > ynumber) {
		nextJ = 1;
	}
	int nextK = k + 1;
	if(nextK > znumber) {
		nextK = 1;
	}
	double ErightY = EfieldZ[i][nextJ][k];
	double EleftY = EfieldZ[i][j][k];
	double ErightZ = EfieldY[i][j][nextK];
	double EleftZ = EfieldY[i][j][k];

	double rotEx = (ErightY - EleftY)/deltaY - (ErightZ - EleftZ)/deltaZ;

	BfieldX[i][j][k] -= speed_of_light_normalized*dt*rotEx;

	alertNaNOrInfinity(BfieldX[i][j][k], "BfieldX = NaN\n");
}

void Simulation::updateBfieldY(int i, int j, int k, double dt) {
	int nextJ = j + 1;
	if(nextJ > ynumber) {
		nextJ = 1;
	}
	int nextK = k + 1;
	if(nextK > znumber) {
		nextK = 1;
	}
	double ErightX = EfieldZ[i+1][j][k];
	double EleftX = EfieldZ[i][j][k];
	double ErightZ = EfieldX[i][j][nextK];
	double EleftZ = EfieldX[i][j][k];

	double rotEy = - (ErightX - EleftX)/deltaX + (ErightZ - EleftZ)/deltaZ;

	BfieldY[i][j][k] -= speed_of_light_normalized*dt*rotEy;

	alertNaNOrInfinity(BfieldY[i][j][k], "BfieldY = NaN\n");
}

void Simulation::updateBfieldZ(int i, int j, int k, double dt) {
	int nextJ = j + 1;
	if(nextJ > ynumber) {
		nextJ = 1;
	}
	int nextK = k + 1;
	if(nextK > znumber) {
		nextK = 1;
	}
	double ErightY = EfieldX[i][nextJ][k];
	double EleftY = EfieldX[i][j][k];
	double ErightX = EfieldY[i+1][j][k];
	double EleftX = EfieldY[i][j][k];

	double rotEz = - (ErightY - EleftY)/deltaY + (ErightX - EleftX)/deltaX;

	BfieldZ[i][j][k] -= speed_of_light_normalized*dt*rotEz;

	alertNaNOrInfinity(BfieldZ[i][j][k], "BfieldZ = NaN\n");
}

void Simulation::updateBoundaries() {
	if(boundaryConditionType == PERIODIC){
		for(int j = 0; j < ynumber; ++j) {
			for(int k = 0; k < znumber; ++k) {
				EfieldY[xnumber][j][k] = EfieldY[0][j][k];
				EfieldZ[xnumber][j][k] = EfieldZ[0][j][k];

				BfieldX[xnumber][j][k] = BfieldX[0][j][k];
			}
		}
	}

	for(int i = 0; i < xnumber + 1; ++i) {
		for(int j = 0; j < ynumber + 1; ++j) {
			if(i < xnumber){
				EfieldX[i][j][znumber] = EfieldX[i][j][0];
				if(j < ynumber) {
					BfieldZ[i][j][znumber] = BfieldZ[i][j][0];
				}
			}

			if(j < ynumber){
				EfieldY[i][j][znumber] = EfieldY[i][j][0];
			}
		}

		for(int k = 0; k < znumber + 1; ++k) {
			if(i < xnumber) {
				EfieldX[i][ynumber][k] = EfieldX[i][0][k];
				if(k < znumber) {
					BfieldY[i][ynumber][k] = BfieldY[i][0][k];
				}
			}
			if(k < znumber){
				EfieldZ[i][ynumber][k] = EfieldZ[i][0][k];
			}
		}
	}
}

double Simulation::evaluateDivE(int i, int j, int k) {
	int prevJ = j - 1;
	if(prevJ < 0) {
		prevJ = ynumber - 1;
	}
	int middleJ = j;
	if(middleJ >= ynumber) {
		middleJ = 0;
	}

	int prevK = k - 1;
	if(prevK < 0) {
		prevK = znumber - 1;
	}
	int middleK = k;
	if(middleK >= znumber) {
		middleK = 0;
	}

	if(boundaryConditionType == SUPERCONDUCTERLEFT){
		if(i == 0) {
			return 0;
		}

		if(i == xnumber) {
			return ((E0.x - EfieldX[i-1][j][k])/deltaX) + ((EfieldY[i][middleJ][k] - EfieldY[i][prevJ][k])/deltaY) + ((EfieldZ[i][j][middleK] - EfieldZ[i][j][prevK])/deltaZ);
		}
	} else if(boundaryConditionType == PERIODIC){
		if(i == 0) {
			return ((EfieldX[i][j][k] - EfieldX[xnumber - 1][j][k])/deltaX) + ((EfieldY[i][middleJ][k] - EfieldY[i][prevJ][k])/deltaY) + ((EfieldZ[i][j][middleK] - EfieldZ[i][j][prevK])/deltaZ);
		} 

		if(i == xnumber) {
			return ((EfieldX[0][j][k] - EfieldX[i-1][j][k])/deltaX) + ((EfieldY[i][middleJ][k] - EfieldY[i][prevJ][k])/deltaY) + ((EfieldZ[i][j][middleK] - EfieldZ[i][j][prevK])/deltaZ);
		}
	}

	return ((EfieldX[i][j][k] - EfieldX[i-1][j][k])/deltaX) + ((EfieldY[i][middleJ][k] - EfieldY[i][prevJ][k])/deltaY) + ((EfieldZ[i][j][middleK] - EfieldZ[i][j][prevK])/deltaZ);
}