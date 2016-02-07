#include <stdio.h>
#include <iostream>
#include <string>
#include <set>

#include "tour.hh"

using namespace std;

int main(int argc, char *argv[]) {
    string inputFilePath = "ConcertScheduling.sampledata";

	int i;
    for (; i < argc; i++) {
    	if (i == 1) {
    		inputFilePath = string(argv[i]);
    		break;
    	}
    } 

    // cout << inputFilePath << endl;

    Tour *tour = Tour::loadCities(inputFilePath);
    Router *router = Router::christofidesRouter();

    vector<City *> cities = tour->route(router);

    if (cities.size() < 1) {
    	cout << "0 km" << endl;

	    delete router;
	    delete tour;
    	return(0);
    }

    City *lastCity = cities[0];
    cout << lastCity->name() << endl;

    double totalDistance = 0;
    set<City *> uniqueCities;
    uniqueCities.insert(lastCity);
    for (i = 1; i < cities.size(); i++) {
    	City *currentCity = cities[i];
    	double distance = lastCity->distance(currentCity);
    	totalDistance += distance;

    	lastCity = currentCity;
        uniqueCities.insert(lastCity);
    	cout << lastCity->name() << endl;
    }

    int intDistance = (totalDistance + 0.5);

    cout << intDistance << " km" << endl;

    delete router;
    delete tour;
    return(0);
}
