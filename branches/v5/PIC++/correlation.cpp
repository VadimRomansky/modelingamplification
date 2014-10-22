#include "stdio.h"
#include "math.h"

#include "util.h"
#include "matrix3d.h"
#include "vector3d.h"
#include "simulation.h"

Vector3d Simulation::correlationEfield(Particle* particle) {
	double x = correlationEfieldX(particle);
	double y = correlationEfieldY(particle);
	double z = correlationEfieldZ(particle);

	return Vector3d(x, y, z);
}

Vector3d Simulation::correlationBfield(Particle* particle) {
	double x = correlationBfieldX(particle);
	double y = correlationBfieldY(particle);
	double z = correlationBfieldZ(particle);

	return Vector3d(x, y, z);
}

double Simulation::correlationEfieldX(Particle* particle) {
	int xcount = truncate(particle->coordinates.x/deltaX);
	int ycount = truncate((particle->coordinates.y/deltaY) + 0.5);
	int zcount = truncate((particle->coordinates.z/deltaZ) + 0.5);

	double Ex = 0;

	for(int i = xcount - 1; i <= xcount + 1; ++i) {
		for(int j = ycount - 1; j <= ycount + 1; ++j) {
			for(int k = zcount - 1; k <= zcount + 1; ++k) {
				Ex += correlationWithExBin(i, j, k, particle)*getEx(i, j, k);
			}
		}
	}

	return Ex;
}

double Simulation::correlationEfieldY(Particle* particle) {
	int xcount = truncate((particle->coordinates.x/deltaX) + 0.5);
	int ycount = truncate(particle->coordinates.y/deltaY);
	int zcount = truncate((particle->coordinates.z/deltaZ) + 0.5);

	double Ey = 0;

	for(int i = xcount - 1; i <= xcount + 1; ++i) {
		for(int j = ycount - 1; j <= ycount + 1; ++j) {
			for(int k = zcount - 1; k <= zcount + 1; ++k) {
				Ey += correlationWithEyBin(i, j, k, particle)*getEy(i, j, k);
			}
		}
	}

	return Ey;
}

double Simulation::correlationEfieldZ(Particle* particle) {
	int xcount = truncate((particle->coordinates.x/deltaX) + 0.5);
	int ycount = truncate((particle->coordinates.y/deltaY) + 0.5);
	int zcount = truncate(particle->coordinates.z/deltaZ);

	double Ez = 0;

	for(int i = xcount - 1; i <= xcount + 1; ++i) {
		for(int j = ycount - 1; j <= ycount + 1; ++j) {
			for(int k = zcount - 1; k <= zcount + 1; ++k) {
				Ez += correlationWithEzBin(i, j, k, particle)*getEz(i, j, k);
			}
		}
	}

	return Ez;
}

double Simulation::correlationBfieldX(Particle* particle) {
	int xcount = truncate((particle->coordinates.x/deltaX) + 0.5);
	int ycount = truncate(particle->coordinates.y/deltaY);
	int zcount = truncate(particle->coordinates.z/deltaZ);

	double Bx = 0;

	for(int i = xcount - 1; i <= xcount + 1; ++i) {
		for(int j = ycount - 1; j <= ycount + 1; ++j) {
			for(int k = zcount - 1; k <= zcount + 1; ++k) {
				Bx += correlationWithBxBin(i, j, k, particle)*getBx(i, j, k);
			}
		}
	}

	return Bx;
}

double Simulation::correlationBfieldY(Particle* particle) {
	int xcount = truncate(particle->coordinates.x/deltaX);
	int ycount = truncate((particle->coordinates.y/deltaY) + 0.5);
	int zcount = truncate(particle->coordinates.z/deltaZ);

	double By = 0;

	for(int i = xcount - 1; i <= xcount + 1; ++i) {
		for(int j = ycount - 1; j <= ycount + 1; ++j) {
			for(int k = zcount - 1; k <= zcount + 1; ++k) {
				By += correlationWithByBin(i, j, k, particle)*getBy(i, j, k);
			}
		}
	}

	return By;
}

double Simulation::correlationBfieldZ(Particle* particle) {
	int xcount = truncate(particle->coordinates.x/deltaX);
	int ycount = truncate(particle->coordinates.y/deltaY);
	int zcount = truncate((particle->coordinates.z/deltaZ) + 0.5);

	double Bz = 0;

	for(int i = xcount - 1; i <= xcount + 1; ++i) {
		for(int j = ycount - 1; j <= ycount + 1; ++j) {
			for(int k = zcount - 1; k <= zcount + 1; ++k) {
				Bz += correlationWithBzBin(i, j, k, particle)*getBz(i, j, k);
			}
		}
	}

	return Bz;
}

double Simulation::correlationWithExBin(int i, int j, int k, Particle* particle) {
	double correlation = 1;

	correlation *= correlationX(i, particle);
	correlation *= correlationShiftY(j, particle);
	correlation *= correlationShiftZ(k, particle);

	return correlation;
}

double Simulation::correlationWithEyBin(int i, int j, int k, Particle* particle) {
	double correlation = 1;

	correlation *= correlationShiftX(i, particle);
	correlation *= correlationY(j, particle);
	correlation *= correlationShiftZ(k, particle);

	return correlation;
}

double Simulation::correlationWithEzBin(int i, int j, int k, Particle* particle) {
	double correlation = 1;

	correlation *= correlationShiftX(i, particle);
	correlation *= correlationShiftY(j, particle);
	correlation *= correlationZ(k, particle);

	return correlation;
}

double Simulation::correlationWithBxBin(int i, int j, int k, Particle* particle) {
	double correlation = 1;

	correlation *= correlationShiftX(i, particle);
	correlation *= correlationY(j, particle);
	correlation *= correlationZ(k, particle);

	return correlation;
}

double Simulation::correlationWithByBin(int i, int j, int k, Particle* particle) {
	double correlation = 1;

	correlation *= correlationX(i, particle);
	correlation *= correlationShiftY(j, particle);
	correlation *= correlationZ(k, particle);

	return correlation;
}

double Simulation::correlationWithBzBin(int i, int j, int k, Particle* particle) {
	double correlation = 1;

	correlation *= correlationX(i, particle);
	correlation *= correlationY(j, particle);
	correlation *= correlationShiftZ(k, particle);

	return correlation;
}

double Simulation::correlationWithGridBin(int i, int j, int k, Particle* particle) {
	double correlation = 1;

	correlation *= correlationX(i, particle);
	correlation *= correlationY(j, particle);
	correlation *= correlationZ(k, particle);

	return correlation;
}

double Simulation::correlationWithShiftGridBin(int i, int j, int k, Particle* particle) {
	double correlation = 1;

	correlation *= correlationShiftX(i, particle);
	correlation *= correlationShiftY(j, particle);
	correlation *= correlationShiftZ(k, particle);

	return correlation;
}

double Simulation::correlationX(int i, Particle* particle) {
	double rightX;
	double leftX;

	if(i < 0) {
		rightX = 0;
		leftX = - deltaX;
	} else if(i >= xnumber){
		rightX = xgrid[xnumber] + deltaX;
		leftX = xgrid[xnumber];
	} else {
		rightX = xgrid[i+1];
		leftX = xgrid[i];
	}

	double correlation = correlationBspline(particle->coordinates.x, particle->dx, leftX, rightX);

	return correlation;
}

double Simulation::correlationShiftX(int i, Particle* particle) {
	double rightX;
	double leftX;

	if(i < 0) {
		rightX = -deltaX/2;
		leftX = - 3*deltaX/2;
	} else if(i > xnumber){
		rightX = xgrid[xnumber] + 3*deltaX/2;
		leftX = xgrid[xnumber] + deltaX/2;
	} else {
		rightX = xgrid[i] + deltaX/2;
		leftX = xgrid[i] - deltaX/2;
	}

	double correlation = correlationBspline(particle->coordinates.x, particle->dx, leftX, rightX);

	return correlation;
}

double Simulation::correlationY(int j, Particle* particle) {
	double rightY;
	double leftY;

	if(j < 0) {
		rightY = ygrid[0];
		leftY = ygrid[0] - deltaY;
	} else if(j >= ynumber) {
		rightY = ygrid[ynumber] + deltaY;
		leftY = ygrid[ynumber];
	} else {
		rightY = ygrid[j + 1];
		leftY = ygrid[j];
	}

	double correlation = correlationBspline(particle->coordinates.y, particle->dy, leftY, rightY);

	return correlation;
}

double Simulation::correlationShiftY(int j, Particle* particle) {
	double rightY;
	double leftY;

	if(j < 0) {
		rightY = ygrid[0] - deltaY/2;
		leftY = ygrid[0] - 3*deltaY/2;
	} else if(j > ynumber) {
		rightY = ygrid[ynumber] + 3*deltaY/2;
		leftY = ygrid[ynumber] + deltaY/2;
	} else {
		rightY = ygrid[j] + deltaY/2;
		leftY = ygrid[j] - deltaY/2;
	}

	double correlation = correlationBspline(particle->coordinates.y, particle->dy, leftY, rightY);

	return correlation;
}

double Simulation::correlationZ(int k, Particle* particle) {
	double rightZ;
	double leftZ;

	if(k < 0) {
		rightZ = zgrid[0];
		leftZ = zgrid[0] - deltaZ;
	} else if(k >= znumber) {
		rightZ = zgrid[znumber] + deltaZ;
		leftZ = zgrid[znumber];
	} else {
		rightZ = zgrid[k + 1];
		leftZ = zgrid[k];
	}

	double correlation = correlationBspline(particle->coordinates.z, particle->dz, leftZ, rightZ);

	return correlation;
}

double Simulation::correlationShiftZ(int k, Particle* particle) {
	double rightZ;
	double leftZ;

	if(k < 0) {
		rightZ = zgrid[0] - deltaZ/2;
		leftZ = zgrid[0] - 3*deltaZ/2;
	} else if(k > znumber) {
		rightZ = zgrid[znumber] + 3*deltaZ/2;
		leftZ = zgrid[znumber] + deltaZ/2;
	} else {
		rightZ = zgrid[k] + deltaZ/2;
		leftZ = zgrid[k] - deltaZ/2;
	}

	double correlation = correlationBspline(particle->coordinates.z, particle->dz, leftZ, rightZ);

	return correlation;
}

double Simulation::correlationBspline(double const& x, double const& dx, double const& leftx, double const& rightx) {
	if (rightx < leftx) {
		printf("rightx < leftx\n");
		exit(0);
	}
	if (dx > rightx - leftx) {
		printf("dx < rightx - leftx\n");
		exit(0);
	}

	double correlation = 0;

	if (x < leftx - dx)
		return 0;
	if (x > rightx + dx)
		return 0;

	if (x < leftx - dx/2) {
		correlation = 2*cube(x + dx - leftx)/(3*cube(dx));
	} else if(x < leftx){
		correlation = (1.0/12.0) + ((x + dx/2 - leftx)/dx) - 2*(cube(dx/2) - cube(leftx - x))/(3*cube(dx));
	} else if (x > rightx + dx/2) {
		correlation = 2*cube(rightx - (x - dx))/(3*cube(dx));
	} else if(x > rightx){
		correlation = (1.0/12.0) + ((-(x - dx/2) + rightx)/dx) - 2*(cube(dx/2) - cube(x - rightx))/(3*cube(dx));
	} else if (x < leftx + dx/2) {
		correlation = 0.5 + ((x - leftx)/dx) - 2*(cube(x - leftx))/(3*cube(dx));
	} else if(x < leftx + dx){
		correlation = 11.0/12.0 + 2*(cube(dx/2) - cube(leftx - (x - dx)))/(3*cube(dx));
	} else if (x > rightx - dx/2) {
		correlation = 0.5 + ((rightx - x)/dx) - 2*(cube(rightx - x))/(3*cube(dx));
	} else if(x > rightx - dx) {
		correlation = 11.0/12.0 + 2*(cube(dx/2) - cube(x + dx - rightx))/(3*cube(dx));
	}else {
		correlation = 1;
	}

	return correlation;
}

double Simulation::getEx(int i, int j, int k) {
	if( i < 0) {
		i = 0;
	}
	if(i >= xnumber) {
		return E0.x;
	}

	if(j < 0) {
		j = ynumber - 1;
	}
	if(j > ynumber) {
		j = 1;
	}

	if(k < 0) {
		k = znumber - 1;
	}
	if(k > znumber) {
		k = 1;
	}

	return EfieldX[i][j][k];
}

double Simulation::getEy(int i, int j, int k) {
	if( i < 0) {
		i = 0;
	}
	if(i > xnumber) {
		return E0.y;
	}

	if(j < 0) {
		j = ynumber - 1;
	}
	if(j >= ynumber) {
		j = 0;
	}

	if(k < 0) {
		k = znumber - 1;
	}
	if(k > znumber) {
		k = 1;
	}

	return EfieldY[i][j][k];
}

double Simulation::getEz(int i, int j, int k) {
	if( i < 0) {
		i = 0;
	}
	if(i > xnumber) {
		return E0.z;
	}

	if(j < 0) {
		j = ynumber - 1;
	}
	if(j > ynumber) {
		j = 1;
	}

	if(k < 0) {
		k = znumber - 1;
	}
	if(k >= znumber) {
		k = 0;
	}

	return EfieldZ[i][j][k];
}

double Simulation::getBx(int i, int j, int k) {
	if( i < 0) {
		i = 0;
	}
	if(i > xnumber) {
		return B0.x;
	}

	if(j < 0) {
		j = ynumber - 1;
	}
	if(j >= ynumber) {
		j = 0;
	}

	if(k < 0) {
		k = znumber - 1;
	}
	if(k >= znumber) {
		k = 0;
	}

	return BfieldX[i][j][k];
}

double Simulation::getBy(int i, int j, int k) {
	if( i < 0) {
		i = 0;
	}
	if(i >= xnumber) {
		return B0.y;
	}

	if(j < 0) {
		j = ynumber - 1;
	}
	if(j > ynumber) {
		j = 1;
	}

	if(k < 0) {
		k = znumber - 1;
	}
	if(k >= znumber) {
		k = 0;
	}

	return BfieldY[i][j][k];
}

double Simulation::getBz(int i, int j, int k) {
	if( i < 0) {
		i = 0;
	}
	if(i >= xnumber) {
		return B0.z;
	}

	if(j < 0) {
		j = ynumber - 1;
	}
	if(j >= ynumber) {
		j = 0;
	}

	if(k < 0) {
		k = znumber - 1;
	}
	if(k > znumber) {
		k = 1;
	}

	return BfieldZ[i][j][k];
}