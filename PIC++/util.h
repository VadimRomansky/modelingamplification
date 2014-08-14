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

double coordinateDifference(double* const a, double* const b, double dt);
void solveSpecialMatrix(double** const leftHalf, double* const rightPart, double* const output);

#endif
