// City.hh
#ifndef CITY_H
#define CITY_H

#include <string>

using namespace std;

class City {
	string _name;
	double _latitude;
	double _longitude;

public:
	string name() const { return _name; }
	double latitude() const { return _latitude; }
	double longitude() const { return _longitude; }
	double distance(City *other);

	City(string name, double latitude, double longitude);
	~City();
};

#endif