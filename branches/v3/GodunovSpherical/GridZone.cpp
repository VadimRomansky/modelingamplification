#include <list>
#include "constants.h"
#include "GridZone.h"

GridZone::GridZone(double central, double left, double right, int t){
	centralPoint = central;
	leftBound = left;
	rightBound = right;
	type = t;
}

GridZone::~GridZone(){
	leftPoints.clear();
	rightPoints.clear();
}

void GridZone::addPoint(){
	double r;
	if(leftPoints.size() == 0){
		r = gridExpLevel*centralPoint + (1 - gridExpLevel)*leftBound;
		leftPoints.push_front(r);
		return;
	}
	if(rightPoints.size() == 0){
		r = gridExpLevel*centralPoint + (1 - gridExpLevel)*rightBound;
		rightPoints.push_front(r);
		return;
	}
	if(leftPoints.size() <= rightPoints.size()){
		r = gridExpLevel*centralPoint + (1 - gridExpLevel)*leftPoints.front();
		leftPoints.push_front(r);
	} else {
		r = gridExpLevel*centralPoint + (1 - gridExpLevel)*rightPoints.front();
		rightPoints.push_front(r);
	}
}