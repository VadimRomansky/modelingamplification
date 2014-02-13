#ifndef _GRID_ZONE_H_
#define _GRID_ZONE_H_

#include <list>

class GridZone{

public:

	double centralPoint;
	double leftBound;
	double rightBound;

	int type;

	std::list<double> leftPoints;
	std::list<double> rightPoints;


	GridZone(double left, double right, int t);
	~GridZone();
	void addPoint();
	void addPoint(double x);
	void addPointLeft(double x);
	void addPointRight(double x);
	void addPoints(int count);
	//void addPoints(int leftCount, int rightCount);

};

#endif