#include "stdafx.h"
#include "BinFlux.h"

BinFlux::BinFlux(){
		fluxR1 = 0;
		fluxR2 = 0;
}
BinFlux::BinFlux(double fr1,double fr2){
		fluxR1 = fr1;
		fluxR2 = fr2;
}
void BinFlux::reset(){
	fluxR1 = 0;
	fluxR2 = 0;
}
