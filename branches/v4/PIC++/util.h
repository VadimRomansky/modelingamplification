#ifndef UTIL_H
#define UTIL_H

double power(double v, double p);
double sqr(double v);
double cube(double v);
double max2(double a, double b);
double min2(double a, double b);
void alertNaNOrInfinity(double value, const char* s);
void alertNotPositive(double value, const char* s);
void alertNegative(double value, const char* s);

void solveSpecialMatrix(double** leftHalf, double* rightPart, double* output);

#endif
