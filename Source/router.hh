// router.hh
#ifndef ROUTER_H
#define ROUTER_H

#include <vector>

#include "city.hh"

using namespace std;

class Router {
public:
	virtual vector<City *> route(vector<City *> cities) {
		return cities;
	}
	static Router * asSpecifiedRouter();
	static Router * nearestNeighbourRouter();
	static Router * christofidesRouter();
};

#endif