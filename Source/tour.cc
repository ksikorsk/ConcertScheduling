#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>

#include "tour.hh"

Tour::Tour(vector<City *> cities) {
	_cities = cities;
}

Tour::~Tour() {
	while (_cities.size() > 0) {
		City *city = _cities[0];
		_cities.erase(_cities.begin());
		delete city;
	}
}

// Split the argument string by the character delimiter
// removing white space and empty entries
vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
    	item.erase(0, item.find_first_not_of(" \n\r\t"));        
    	item.erase(item.find_last_not_of(" \n\r\t")+1);
    	if (!item.empty()) {
	        elems.push_back(item);
    	}
    }
    return elems;
}

// Split the argument string by the character delimiter
// removing white space and empty entries
vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

Tour * Tour::loadCities(string inputPath) {
	ifstream inputFile(inputPath);

	vector<City *> cities;

	if (!inputFile || !inputFile.is_open()) {
		return NULL;
	}

	string line;
	while ( getline(inputFile, line) ) {
	  	vector<string> elems = split(line, ',');
	  	if (elems.size() == 5) {
	  		string name = elems[0];
	  		int latDegrees = stoi(elems[1]);
	  		int latMintutes = stoi(elems[2]);
	  		double latitude = latDegrees + latMintutes / 60.0;

	  		int longDegrees = stoi(elems[3]);
	  		int longMintutes = stoi(elems[4]);
	  		double longitude = longDegrees + longMintutes / 60.0;

	  		City *city = new City(name, latitude, longitude);
	  		cities.push_back(city);
	  	}
    }
	inputFile.close();

    Tour *tour = new Tour(cities);
    return tour;
}