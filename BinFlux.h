#ifndef BINFLUX_H
#define BINFLUX_H

class BinFlux{
public:
	double fluxR1;
	double fluxR2;
	double fluxPhi1;
	double fluxPhi2;
	double fluxTheta1;
	double fluxTheta2;
	BinFlux();
	BinFlux(double fr1,double fr2,double ftheta1,double ftheta2,double fphi1,double fphi2);
	void reset();
};



#endif
