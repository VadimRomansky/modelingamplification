#include "stdafx.h"
#include "BinFlux.h"

BinFlux::BinFlux(){
		fluxR1 = 0;
		fluxR2 = 0;
		fluxPhi1 = 0;
		fluxPhi2 = 0;
		fluxTheta1 = 0;
		fluxTheta2 = 0;
}
BinFlux::BinFlux(double fr1,double fr2,double ftheta1,double ftheta2,double fphi1,double fphi2){
		fluxR1 = fr1;
		fluxR2 = fr2;
		fluxPhi1 = ftheta1;
		fluxPhi2 = ftheta2;
		fluxTheta1 = fphi1;
		fluxTheta2 = fphi2;	
}
void BinFlux::reset(){
	fluxR1 = 0;
	fluxR2 = 0;
	fluxTheta1 = 0;
	fluxTheta2 = 0;
	fluxPhi1 = 0;
	fluxPhi2 = 0;
}
