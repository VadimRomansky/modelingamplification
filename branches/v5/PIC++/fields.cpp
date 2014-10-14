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

	if(i == xnumber) {
		EfieldX[i][j][k] = E0.x;
		return;
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

	EfieldX[i][j][k] += (speed_of_light_normalized*rotBx -4*pi*electricFluxX[i][j][k])*dt;
}

void Simulation::updateEfieldY(int i, int j, int k, double dt) {
	if(i == 0) {
		EfieldY[i][j][k] = 0;
		return;
	}
	if(i == xnumber) {
		EfieldY[i][j][k] = E0.y;
		return;
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

	double BrightX = BfieldZ[i][middleJ][middleK];
	double BleftX = BfieldZ[i-1][middleJ][middleK];

	double BrightZ = BfieldX[i][middleJ][middleK];
	double BleftZ = BfieldX[i][middleJ][prevK];

	double rotBy = -(BrightX - BleftX)/deltaX + (BrightZ - BleftZ)/deltaZ;

	EfieldY[i][j][k] += (speed_of_light_normalized*rotBy - 4*pi*electricFluxY[i][j][k])*dt;
}

void Simulation::updateEfieldZ(int i, int j, int k, double dt) {
	if(i == 0) {
		EfieldZ[i][j][k] = 0;
		return;
	}
	if(i == xnumber) {
		EfieldZ[i][j][k] = E0.z;
		return;
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

	double BrightX = BfieldY[i][middleJ][middleK];
	double BleftX = BfieldY[i-1][middleJ][middleK];

	double BrightY = BfieldX[i][middleJ][middleK];
	double BleftY = BfieldX[i][prevJ][middleK];

	double rotBy = (BrightX - BleftX)/deltaX - (BrightY - BleftY)/deltaY;

	EfieldZ[i][j][k] += (speed_of_light_normalized*rotBy - 4*pi*electricFluxY[i][j][k])*dt;
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
}

void Simulation::updateBoundaries() {
}