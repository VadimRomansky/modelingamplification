#ifndef UTIL_H
#define UTIL_H

double power(double v, double p);
double sqr(double v);
double cube(double v);
double max2(double a, double b);
double max3(double a, double b, double c);
double min2(double a, double b);
double min3(double a, double b, double c);
void alertNaNOrInfinity(double value, const char* s);
void alertNotPositive(double value, const char* s);
void alertNegative(double value, const char* s);
double uniformDistribution();

double coordinateDifference(double* const a, double* const b, double dt);
void solveSpecialMatrix(double** const leftHalf, double* const rightPart, double* const output);

double McDonaldFunction(double x, double index);
double normalDistribution();
double maxwellDistribution(double temperature);
double maxwellJuttnerDistribution(double temperature, double mass);
double solveInverceJuttnerFunction(double x, double theta, double besselK);
double solveInverceJuttnerFunction(double x, double theta, double besselK, double left, double right);
double maxwellJuttnerFunction(double gamma, double theta, double besselK);
double maxwellJuttnerIntegral(double gamma, double theta, double besselK);

#endif
