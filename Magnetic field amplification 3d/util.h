#ifndef UTIL_H
#define UTIL_H
#include "vector3d.h"

double power(double v, double p);
double sqr(double v);
int lowerInt(double v);
vector3d summVelocity(vector3d v, double u);
bool order(double a, double b, double c);
double angleDelta(double phi1, double phi2);

#endif