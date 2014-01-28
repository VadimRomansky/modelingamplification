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
void outputParticlePath(std::list<Particle*>& list,const char* cosmicRayFileName,const char* notCosmicRayFileName);
void outputRadialProfile(SpaceBin** bins, FILE* outProfile, const int rgridNumber);
void outputProfile(SpaceBin*** bins, const int xgridNumber, const int ygridNumber);
void outputShockWave(std::list<double> points, std::list<double> velocity);
void outputAverageVz(double minp, double maxp, double* averageVz, const char* fileName);
void outputParticles(std::vector<Particle*>& particles, const char* fileName);
//void outputParticlesCount(Xbin** bins, const char* fileName);

#endif