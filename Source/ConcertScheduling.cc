#include <stdio.h>
#include <iostream>
#include <string>
#include <set>

#include "tour.hh"

using namespace std;

enum Algorithm {
    none = 0,
    christofides = 1,
    nearestNeighbour = 2,
    asSpecified = 3
};

int main(int argc, char *argv[]) {
    string inputFilePath = "ConcertScheduling.sampledata";
    Algorithm algorithm = none;

    for (int i = 1; i < argc; i++) {
        string arg = string(argv[i]);
        if (arg == "-i") {
            i++;
            inputFilePath = string(argv[i]);
        } else if (arg == "-r") {
            i++;
            string value = string(argv[i]);
            algorithm = (Algorithm)stoi(value);
        }
    } 

    // cout << inputFilePath << endl;

    Tour *tour = Tour::loadCities(inputFilePath);
    Router *router;
    switch(algorithm) {
        case nearestNeighbour: 
            router = Router::nearestNeighbourRouter();
            break;
        case asSpecified: 
            router = Router::asSpecifiedRouter();
            break;
        default:
            router = Router::christofidesRouter();
            break;
    }

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
    for (int i = 1; i < cities.size(); i++) {
    	City *currentCity = cities[i];
    	double distance = lastCity->distance(currentCity);
    	totalDistance += distance;

    	lastCity = currentCity;
    	cout << lastCity->name() << endl;
    }

    int intDistance = (totalDistance + 0.5);

    cout << intDistance << " km" << endl;

    delete router;
    delete tour;
    return(0);
}
