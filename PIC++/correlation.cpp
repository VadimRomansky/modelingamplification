#include "util.h"
#include "simulation.h"

Vector3d Simulation::correlationTempEfield(Particle* particle) {
	return correlationTempEfield(*particle);
}

Vector3d Simulation::correlationBfield(Particle* particle) {
	return correlationBfield(*particle);
}

Vector3d Simulation::correlationEfield(Particle* particle) {
	return correlationEfield(*particle);
}

Vector3d Simulation::correlationTempEfield(Particle& particle) {
	checkParticleInBox(particle);

	int xcount = trunc(particle.coordinates.x / deltaX + 0.5);
	int ycount = trunc(particle.coordinates.y / deltaY + 0.5);
	int zcount = trunc(particle.coordinates.z / deltaZ + 0.5);

	if (xcount < 0) {
		printf("xcount < 0\n");
		//exit(0);
	}

	if (xcount > xnumber) {
		printf("xcount > xnumber\n");
		//exit(0);
	}

	if (ycount < 0) {
		printf("ycount < 0\n");
		//exit(0);
	}

	if (ycount > ynumber) {
		printf("ycount > ynumber\n");
		//exit(0);
	}

	if (zcount < 0) {
		printf("zcount < 0\n");
		//exit(0);
	}

	if (zcount > znumber) {
		printf("zcount > znumber\n");
		//exit(0);
	}

	Vector3d result = Vector3d(0, 0, 0);

	for (int i = xcount - 1; i <= xcount + 1; ++i) {
		for (int j = ycount - 1; j <= ycount + 1; ++j) {
			for (int k = zcount - 1; k <= zcount + 1; ++k) {
				result = result + correlationFieldWithTempEbin(particle, i, j, k);
			}
		}
	}

	return result;
}

Vector3d Simulation::correlationBfield(Particle& particle) {
	checkParticleInBox(particle);

	int xcount = trunc(particle.coordinates.x / deltaX);
	int ycount = trunc(particle.coordinates.y / deltaY);
	int zcount = trunc(particle.coordinates.z / deltaZ);

	if (xcount < 0) {
		printf("xcount < 0\n");
		//exit(0);
	}

	if (xcount >= xnumber) {
		printf("xcount >= xnumber\n");
		//exit(0);
	}

	if (ycount < 0) {
		printf("ycount < 0\n");
		//exit(0);
	}

	if (ycount >= ynumber) {
		printf("ycount >= ynumber\n");
		//exit(0);
	}

	if (zcount < 0) {
		printf("zcount < 0\n");
		//exit(0);
	}

	if (zcount >= znumber) {
		printf("zcount >= znumber\n");
		//exit(0);
	}

	Vector3d result = Vector3d(0, 0, 0);

	for (int i = xcount - 1; i <= xcount + 1; ++i) {
		for (int j = ycount - 1; j <= ycount + 1; ++j) {
			for (int k = zcount - 1; k <= zcount + 1; ++k) {
				result = result + correlationFieldWithBbin(particle, i, j, k);
			}
		}
	}

	return result;
}

Vector3d Simulation::correlationEfield(Particle& particle) {
	checkParticleInBox(particle);

	int xcount = trunc(particle.coordinates.x / deltaX + 0.5);
	int ycount = trunc(particle.coordinates.y / deltaY + 0.5);
	int zcount = trunc(particle.coordinates.z / deltaZ + 0.5);

	if (xcount < 0) {
		printf("xcount < 0\n");
		//exit(0);
	}

	if (xcount > xnumber) {
		printf("xcount > xnumber\n");
		//exit(0);
	}

	if (ycount < 0) {
		printf("ycount < 0\n");
		//exit(0);
	}

	if (ycount > ynumber) {
		printf("ycount > ynumber\n");
		//exit(0);
	}

	if (zcount < 0) {
		printf("zcount < 0\n");
		//exit(0);
	}

	if (zcount > znumber) {
		printf("zcount > znumber\n");
		//exit(0);
	}

	Vector3d result = Vector3d(0, 0, 0);

	for (int i = xcount - 1; i <= xcount + 1; ++i) {
		for (int j = ycount - 1; j <= ycount + 1; ++j) {
			for (int k = zcount - 1; k <= zcount + 1; ++k) {
				result = result + correlationFieldWithEbin(particle, i, j, k);
			}
		}
	}

	return result;
}


Vector3d Simulation::correlationFieldWithBbin(Particle& particle, int i, int j, int k) {

	Vector3d _field;

	if (i < 0) {
		//_field = Vector3d(0,0,0);
		//note: not zero because reflection
		if (j < 0) {
			if (k < 0) {
				_field = Bfield[0][ynumber - 1][znumber - 1];
			} else if (k >= znumber) {
				_field = Bfield[0][ynumber - 1][0];
			} else {
				_field = Bfield[0][ynumber - 1][k];
			}
		} else if (j >= ynumber) {
			if (k < 0) {
				_field = Bfield[0][0][znumber - 1];
			} else if (k >= znumber) {
				_field = Bfield[0][0][0];
			} else {
				_field = Bfield[0][0][k];
			}
		} else {
			if (k < 0) {
				_field = Bfield[0][j][znumber - 1];
			} else if (k >= znumber) {
				_field = Bfield[0][j][0];
			} else {
				_field = Bfield[0][j][k];
			}
		}
	} else if (i >= xnumber) {
		_field = B0;
	} else {
		if (j < 0) {
			if (k < 0) {
				_field = Bfield[i][ynumber - 1][znumber - 1];
			} else if (k >= znumber) {
				_field = Bfield[i][ynumber - 1][0];
			} else {
				_field = Bfield[i][ynumber - 1][k];
			}
		} else if (j >= ynumber) {
			if (k < 0) {
				_field = Bfield[i][0][znumber - 1];
			} else if (k >= znumber) {
				_field = Bfield[i][0][0];
			} else {
				_field = Bfield[i][0][k];
			}
		} else {
			if (k < 0) {
				_field = Bfield[i][j][znumber - 1];
			} else if (k >= znumber) {
				_field = Bfield[i][j][0];
			} else {
				_field = Bfield[i][j][k];
			}
		}
	}

	double correlation = correlationWithBbin(particle, i, j, k);

	return _field * correlation;
}

Vector3d Simulation::correlationFieldWithEbin(Particle& particle, int i, int j, int k) {

	Vector3d _field = getEfield(i, j, k);

	double correlation = correlationWithEbin(particle, i, j, k);

	return _field * correlation;
}

Vector3d Simulation::correlationFieldWithTempEbin(Particle& particle, int i, int j, int k) {

	Vector3d _field = getTempEfield(i, j, k);

	double correlation = correlationWithEbin(particle, i, j, k);

	return _field * correlation;
}

double Simulation::correlationWithBbin(Particle& particle, int i, int j, int k) {
	if (! particleCrossBbin(particle, i, j, k))
		return 0.0;

	double x = particle.coordinates.x;
	double y = particle.coordinates.y;
	double z = particle.coordinates.z;

	double leftx;
	double rightx;
	double lefty;
	double righty;
	double leftz;
	double rightz;

	if (i < 0) {
		leftx = particle.coordinates.x - 2 * deltaX;
		rightx = xgrid[0];
	} else if (i >= xnumber) {
		leftx = xgrid[xnumber];
		rightx = particle.coordinates.x + 2 * deltaX;
	} else {
		leftx = xgrid[i];
		rightx = xgrid[i + 1];

	}


	if (j < 0) {
		lefty = ygrid[0] - deltaY;
		righty = ygrid[0];
	} else if (j >= ynumber) {
		lefty = ygrid[ynumber];
		righty = ygrid[ynumber] + deltaY;
	} else if (j == 0) {
		if (particle.coordinates.y - particle.dy > ygrid[1]) {
			lefty = ygrid[ynumber];
			righty = ygrid[ynumber] + deltaY;
		} else {
			lefty = ygrid[0];
			righty = ygrid[1];
		}
	} else if (j == ynumber - 1) {
		if (particle.coordinates.y + particle.dy < ygrid[ynumber - 1]) {
			lefty = ygrid[0] - deltaY;
			righty = ygrid[0];
		} else {
			lefty = ygrid[ynumber - 1];
			righty = ygrid[ynumber];
		}
	} else {
		lefty = ygrid[j];
		righty = ygrid[j + 1];
	}

	if (k < 0) {
		leftz = zgrid[0] - deltaZ;
		rightz = zgrid[0];
	} else if (k >= znumber) {
		leftz = zgrid[znumber];
		rightz = zgrid[znumber] + deltaZ;
	} else if (k == 0) {
		if (particle.coordinates.z - particle.dz > zgrid[1]) {
			leftz = zgrid[ynumber];
			rightz = zgrid[ynumber] + deltaZ;
		} else {
			leftz = zgrid[0];
			rightz = zgrid[1];
		}
	} else if (k == znumber - 1) {
		if (particle.coordinates.z + particle.dz < zgrid[znumber - 1]) {
			leftz = zgrid[0] - deltaZ;
			rightz = zgrid[0];
		} else {
			leftz = zgrid[znumber - 1];
			rightz = zgrid[znumber];
		}
	} else {
		leftz = zgrid[k];
		rightz = zgrid[k + 1];
	}


	double correlation = 1;

	double correlationx = correlationBspline(x, particle.dx, leftx, rightx);
	double correlationy = correlationBspline(y, particle.dy, lefty, righty);
	double correlationz = correlationBspline(z, particle.dz, leftz, rightz);

	correlation = correlationx * correlationy * correlationz;

	return correlation;
}

double Simulation::correlationWithEbin(Particle& particle, int i, int j, int k) {
	double x = particle.coordinates.x;
	double y = particle.coordinates.y;
	double z = particle.coordinates.z;

	double leftx;
	double rightx;
	double lefty;
	double righty;
	double leftz;
	double rightz;

	if (i < 0) {
		leftx = particle.coordinates.x - 2 * deltaX;
		rightx = xgrid[0] - deltaX / 2;
	} else if (i > xnumber) {
		leftx = xgrid[xnumber] + deltaX / 2;
		rightx = particle.coordinates.x + 2 * deltaX;
		/*} else if(i == 0){
		////note: needs because in the middle of 0 Ebin is wall
		leftx = 0;
		rightx = deltaX/2;*/
	} else {
		leftx = xgrid[i] - deltaX / 2;
		rightx = xgrid[i] + deltaX / 2;
	}


	if (j < 0) {
		lefty = particle.coordinates.y - 2 * deltaY;
		righty = ygrid[0] - deltaY / 2;
	} else if (j > ynumber) {
		lefty = ygrid[ynumber] + deltaX / 2;
		righty = particle.coordinates.y + 2 * deltaY;
	} else {
		lefty = ygrid[j] - deltaY / 2;
		righty = ygrid[j] + deltaY / 2;
	}

	if (k < 0) {
		leftz = particle.coordinates.z - 2 * deltaZ;
		rightz = zgrid[0] - deltaZ/2;
	} else if (k > znumber) {
		leftz = zgrid[znumber] + deltaZ / 2;
		rightz = particle.coordinates.z + 2 * deltaZ;
	} else {
		leftz = zgrid[k] - deltaZ / 2;
		rightz = zgrid[k] + deltaZ / 2;
	}


	double correlation = 1;

	double correlationx = correlationBspline(x, particle.dx, leftx, rightx);
	double correlationy = correlationBspline(y, particle.dy, lefty, righty);
	double correlationz = correlationBspline(z, particle.dz, leftz, rightz);

	correlation = correlationx * correlationy * correlationz;

	return correlation;
}

double Simulation::correlationBspline(const double& x, const double& dx, const double& leftx, const double& rightx) {

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

	if (x < leftx) {
		correlation = sqr(x + dx - leftx) / 2;
	} else if (x > rightx) {
		correlation = sqr(rightx - (x - dx)) / 2;
	} else if (x < leftx + dx) {
		correlation = sqr(dx) - sqr(leftx - (x - dx)) / 2;
	} else if (x > rightx - dx) {
		correlation = sqr(dx) - sqr(x + dx - rightx) / 2;
	} else {
		correlation = dx * dx;
	}

	correlation /= dx * dx;

	return correlation;
}