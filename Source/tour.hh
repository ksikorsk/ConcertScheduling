// City.hh
#ifndef TOUR_H
#define TOUR_H

#include <vector>
#include <string>

#include "router.hh"

using namespace std;

class Tour {
	vector<City *> _cities;
	Tour(vector<City *> cities);

public:
	vector<City *> cities() { return _cities; }

	~Tour();
	static Tour * loadCities(string inputPath);
	vector<City *> route(Router * router) {
		return router->route(_cities);
	}
};

#endif