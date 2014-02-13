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
		if(r < leftPoints.front()){
			printf("aaa\n");
		}
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
		/*if(r < leftPoints.front()){
			printf("aaa\n");
		}*/
		leftPoints.push_front(r);
	} else {
		r = x*centralPoint + (1 - x)*rightPoints.front();
		rightPoints.push_front(r);
	}
}

void GridZone::addPointLeft(double x){
	double r;
	if(leftPoints.size() == 0){
		r = x*centralPoint + (1 - x)*leftBound;
		leftPoints.push_front(r);
		return;
	}
	r = x*centralPoint + (1 - x)*leftPoints.front();
	/*if(r < leftPoints.front()){
		printf("aaa\n");
	}*/
	leftPoints.push_front(r);
}

void GridZone::addPointRight(double x){
	double r;
	if(rightPoints.size() == 0){
		r = x*centralPoint + (1 - x)*rightBound;
		rightPoints.push_front(r);
		return;
	}
	r = x*centralPoint + (1 - x)*rightPoints.front();
	rightPoints.push_front(r);
}

void GridZone::addPoints(int count){
	if(count <= 0) return;
	if(type != 0){
		double leftR = centralPoint - leftBound;
		double rightR = rightBound - centralPoint;
		double leftX = gridExpLevel;
		double rightX = gridExpLevel;
		int leftCount = 0;
		int rightCount = 0;
		if(leftR < minDeltaR){
			rightCount = count;
		} else if(rightR < minDeltaR){
			leftCount = count;
		} else {
			rightCount = min(rightR/minDeltaR, rightR*count/(rightBound - leftBound));
			leftCount = count - rightCount;
			leftX = min(1 - exp(log(minDeltaR/leftR)*leftCount), gridExpLevel);
			rightX = min(1 - exp(log(minDeltaR/rightR)*rightCount), gridExpLevel);
		}

		if(leftX < 0) {
			leftX = gridExpLevel;
			printf("deltaR < minDeltaR\n");
		}
		if(rightX < 0) {
			rightX = gridExpLevel;
			printf("deltaR < minDeltaR\n");
		}
		for(int i = 0; i < leftCount; ++i){
			addPointLeft(leftX);
		}
		for(int i = 0; i < rightCount; ++i){
			addPointRight(rightX);
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
				/*if(r < leftPoints.front()){
					printf("aaa\n");
				}*/
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

/*void GridZone::addPoints(int leftCount, int rightCount){
	if((leftCount <= 0) && (rightCount <= 0)) return;
	if(type != 0){
		double R = min(centralPoint - leftBound, rightBound - centralPoint);
		double x = min(1 - exp(log(minDeltaR/R)*2.0/count), gridExpLevel);
		if(x < 0) {
			x = gridExpLevel;
			printf("deltaR < minDeltaR\n");
		}
		for(int i = 0; i < count; ++i){
			addPoint(x);
		}
	} else {
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
}*/