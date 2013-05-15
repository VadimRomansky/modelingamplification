#ifndef OUTPUT_H
#define OUTPUT_H
#include "simulation.h"
#include <list>

void output(Simulation& simulation);
void outputPDF(std::list <Particle*>& list,const char* fileName);
void outputPDF(std::vector <Particle*>& list,const char* fileName);
void outputRPDF(std::list <Particle*>& list,const char* fileName);
void outputEnergyPDF(std::list <Particle*>& list,const char* fileName);
void outputEnergyPDF(std::vector <Particle*>& list,const char* fileName);
void outputStartPDF(std::list <Particle*>& l,const char* fileName, Simulation& simulation,double minp,double maxp);
void outputTurbulenceSpectrum(double* w, const char* fileName,double minK, double maxK);
void outputMagneticField(SpaceBin** bins, const char* fileName);
void outputParticlePath(std::list<Particle*>& list,const char* cosmicRayFileName,const char* notCosmicRayFileName);
void outputRadialProfile(SpaceBin** bins, FILE* outProfile, const int rgridNumber);
void outputShockWave(std::list<double> points, std::list<double> velocity);
void outputAverageVz(double minp, double maxp, double* averageVz, const char* fileName);
void outputParticles(std::vector<Particle*>& particles, const char* fileName);
void outputDistribution(std::vector<Particle*>& partcles, FILE* file, double maxP);
//void outputParticlesCount(Xbin** bins, const char* fileName);

#endif