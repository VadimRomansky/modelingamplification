#include "stdafx.h"
#include <list>
#include "constants.h"
#include "util.h"
#include "GridZone.h"

GridZone::GridZone(double left, double right, int t){
	leftBound = left;
	rightBound = right;
	centralPoint = (left + right)/2;
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

void GridZone::addPoint(double x){
	double r;
	if(leftPoints.size() == 0){
		r = x*centralPoint + (1 - x)*leftBound;
		leftPoints.push_front(r);
		return;
	}
	if(rightPoints.size() == 0){
		r = x*centralPoint + (1 - x)*rightBound;
		rightPoints.push_front(r);
		return;
	}
	if(leftPoints.size() <= rightPoints.size()){
		r = x*centralPoint + (1 - x)*leftPoints.front();
		leftPoints.push_front(r);
	} else {
		r = x*centralPoint + (1 - x)*rightPoints.front();
		rightPoints.push_front(r);
	}
}

void GridZone::addPoint(int count){
	if(count <= 0) return;
	if(type != 0){
		double R = min(centralPoint - leftBound, rightBound - centralPoint);
		double x = min(1 - exp(log(minDeltaR/R)*2/count), gridExpLevel);
		for(int i = 0; i < count; ++i){
			addPoint(x);
		}
	} else {
		int leftCount = count/2;
		int rightCount = count - leftCount;
		double left;
		double right;
		if(leftPoints.size() == 0){
			left = leftBound;
		} else {
			left = leftPoints.front();
		}
		if(rightPoints.size() == 0){
			right  = rightBound;
		} else {
			right = rightPoints.front();
		}
		if(leftCount > 0){
			double dr = (centralPoint - left)/(leftCount + 1);
			double r = left + dr;
			for(int i = 1; i <= leftCount; ++i){
				leftPoints.push_front(r);
				r = r + dr;
			}
		}

		if(rightCount > 0){
			double dr = (right - centralPoint)/(rightCount + 1);
			double r = right - dr;
			for(int i = 1; i <= rightCount; ++i){
				rightPoints.push_front(r);
				r = r - dr;
			}
		}
	}
}