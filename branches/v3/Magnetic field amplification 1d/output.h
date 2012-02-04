#ifndef OUTPUT_H
#define OUTPUT_H
#include "simulation.h"
#include <list>

void output(Simulation& simulation);
void outputPDF(std::list <Particle*> list,const char* fileName);
void outputStartPDF(std::list <Particle*> l,const char* fileName, Simulation& simulation,double minp,double maxp);
void outputTurbulenceSpectrum(double* w, const char* fileName,double minK, double maxK);
void outputMagneticField(SpaceBin**** bins, const char* fileName);
void outputParticlePath(std::list<Particle*> list,const char* cosmicRayFileName,const char* notCosmicRayFileName);
void outputRadialProfile(SpaceBin**** bins, int thetaNumber, int phiNumber, FILE* outProfile, double* averageVelocity);
void outputPhiProfile(SpaceBin**** bins, int rNumber, int thetaNumber);
void outputThetaProfile(SpaceBin**** bins, int rNumber, int phiNumber);
void outputShockWave(std::list<double> points, std::list<double> velocity);
//void outputParticlesCount(Xbin** bins, const char* fileName);

#endif