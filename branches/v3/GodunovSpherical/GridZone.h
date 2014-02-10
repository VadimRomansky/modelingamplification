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


	GridZone(double central, double left, double right, int t);
	~GridZone();
	void addPoint();

};

#endif