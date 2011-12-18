#ifndef BINFLUX_H
#define BINFLUX_H

class BinFlux{
public:
	double fluxR1;
	double fluxR2;
	BinFlux();
	BinFlux(double fr1,double fr2);
	void reset();
};



#endif
