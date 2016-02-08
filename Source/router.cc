#include "router.hh"

#include <iostream>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <stack>

Router * Router::asSpecifiedRouter() {
	return new Router();
}

class NearestNeighbourRouter : public Router {
public:
	vector<City *> route(vector<City *> cities) {
		if (cities.size() < 3) {
			return cities;
		} 

		vector<City *> remainingCities(cities);
		vector<City *> routedCities;
		routedCities.reserve(cities.size());

		City *lastCity = remainingCities[0];
		remainingCities.erase(remainingCities.begin());
		routedCities.push_back(lastCity);

		while(!remainingCities.empty()) {
			int minIndex = numeric_limits<int>::max();
			double minDistance = numeric_limits<double>::max();

			for (int i = 0; i < remainingCities.size(); i++) {
				City *currentCity = remainingCities[i];
				double distance = currentCity->distance(lastCity);

				if (distance < minDistance) {
					minIndex = i;
					minDistance = distance;
				}
			}

			City *minCity = remainingCities[minIndex];
			remainingCities.erase(remainingCities.begin() + minIndex);
			routedCities.push_back(minCity);
		}

		return routedCities;
	}
};

Router * Router::nearestNeighbourRouter() {
	return new NearestNeighbourRouter();
}

struct Edge {
	int _index1, _index2;
	double _weight;

public:
	int index1() { return _index1; }
	int index2() { return _index2; }
	double weight() { return _weight; }

	Edge(int index1, int index2, double weight) {
		_index1 = index1;
		_index2 = index2;
		_weight = weight;
	}

	bool hasIndex(int index) {
		return _index1 == index || _index2 == index;
	}

	bool hasIndex(Edge * edge) {
		return hasIndex(edge->index1()) || hasIndex(edge->index2());
	}

	~Edge() {
	}
};

bool edgeCompare(Edge * a, Edge * b) {
	return a->weight() < b->weight();
}

// References:
// Krushkal's Algorithm
// http://www.geeksforgeeks.org/greedy-algorithms-set-2-kruskals-minimum-spanning-tree-mst/
// Github: beckysag, Traveling Salesperson Problem
// https://github.com/beckysag/traveling-salesman
class ChristofidesRouter : public Router {
	// A structure to represent a subset for union-find
	struct Subset {
	public:
		int parent;
	    int rank;

	    Subset(int parent) {
	    	this->parent = parent;
	    	this->rank = 0;
	    }

	    ~Subset() {
	    }
	};

	// A utility function to find set of an element i
	// (uses path compression technique)
	int findSubset(vector<Subset *> subsets, int index) {
	    // find root and make root as parent of i (path compression)
	    if (subsets[index]->parent != index) {
	        subsets[index]->parent = findSubset(subsets, subsets[index]->parent);
	    }
	 
	    return subsets[index]->parent;
	}

	// A function that does union of two sets of x and y
	// (uses union by rank)
	void unionSubset(vector<Subset *> subsets, int x, int y) {
	    int xroot = findSubset(subsets, x);
	    int yroot = findSubset(subsets, y);
	 
	    // Attach smaller rank tree under root of high rank tree
	    // (Union by Rank)
	    if (subsets[xroot]->rank < subsets[yroot]->rank) {
	        subsets[xroot]->parent = yroot;
	    }
	    else if (subsets[xroot]->rank > subsets[yroot]->rank) {
	        subsets[yroot]->parent = xroot;
	    }
	 
	    // If ranks are same, then make one as root and increment
	    // its rank by one
	    else
	    {
	        subsets[yroot]->parent = xroot;
	        subsets[xroot]->rank++;
	    }
	}

	// void printEdge(vector<City *> cities, Edge * edge) {
	// 	City * city1 = cities[edge->index1()];
	// 	City * city2 = cities[edge->index2()];

	// 	cout << city1->name() << " to " << city2->name() << endl;
	// }

	// 1) Create a minimum spanning tree T of G.
	//    Kruskal's Algorithm
	vector<Edge *> kruskalsMst(vector<City *> vertices, vector<Edge *> edges) {
		// cout << "Step 1: Kruskal's MST" << endl;

		vector<Subset *> subsets;
		for(int i = 0; i < vertices.size(); i++) {
			subsets.push_back(new Subset(i));
		}

		// cout << "Step 1a: Sorting edges" << endl;

		// Step 1:  Sort all the edges in non-decreasing order of their weight
		// If we are not allowed to change the given graph, 
		// we can create a copy of array of edges
		sort(edges.begin(), edges.end(), edgeCompare);

		// cout << "Step 1b: Processing edges " << vertices.size() << endl;

		// Number of edges to be taken is equal to V-1
		int index = 0;
		vector<Edge *> mstEdges;
		while (mstEdges.size() < (vertices.size() - 1)) {
			// Step 2: Pick the smallest edge. And increment the index
			// for next iteration

			Edge * edge = edges[index++];

			int x = findSubset(subsets, edge->index1());
			int y = findSubset(subsets, edge->index2());

			// // If including this edge does't cause cycle, include it
			// // in result and increment the index of result for next edge

			if (x != y) {
				mstEdges.push_back(edge);
				unionSubset(subsets, x, y);

				// printEdge(vertices, edge);
				continue;
			}
		}

		while (!subsets.empty()) {
			Subset * subset = subsets[0];
			delete subset;
			subsets.erase(subsets.begin());
		}

		return mstEdges;		
	}

	// 2) Let O be the set of vertices with odd degree in T. 
	//    By the handshaking lemma, O has an even number of vertices.
	// 3) Find a minimum-weight perfect matching M in the induced
	//    subgraph given by the vertices from O.
	vector<Edge *> perfectMatching(vector<City *> vertices, vector<Edge *> edges, vector<Edge *> mstEdges) {
		// cout << "Step 2: Set of verticies O with odd degree" << endl;

		unordered_map<int, int> countHash(vertices.size());
		for (int i = 0; i < mstEdges.size(); i++) {
			Edge * edge = mstEdges[i];
			countHash[edge->index1()]++;
			countHash[edge->index2()]++;
		}

		// cout << "Step 3: minimum-weight perfect matching M from O" << endl;

		// cout << "Step 3a: Form the subgraph of G using only the vertices of O" << endl;

		vector<Edge *> subEdges;
		for (int i = 0; i < edges.size(); i++) {
			Edge * edge = edges[i];
			bool city1Odd = (countHash[edge->index1()] % 2) > 0;
			bool city2Odd = (countHash[edge->index2()] % 2) > 0;
			if (city1Odd && city2Odd) {
				subEdges.push_back(edge);
				// printEdge(vertices, edge);
			}
		}

		// cout << "Step 3b: Construct a perfect matching M in this subgraph using greedy (not min) algorithm" << endl;

		// for each odd node
		vector<Edge *> perfectMatching;
		while (!subEdges.empty()) {
			double minWeight = numeric_limits<double>::max();
			Edge * minEdge = NULL;
			for (int i = 0; i < subEdges.size(); i++) {
				Edge * edge = subEdges[i];
				if (edge->weight() < minWeight) {
					minWeight = edge->weight();
					minEdge = edge;
				}
			}

			perfectMatching.push_back(minEdge);

			// Remove all the edges connected to minEdge (including minEdge)
			for (int i = 0; i < subEdges.size(); i++) {
				Edge * edge = subEdges[i];
				if (edge->hasIndex(minEdge)) {
					subEdges.erase(subEdges.begin() + i);
					i--;
				}
			}
		}

		for (int i = 0; i < perfectMatching.size(); i++) {
			Edge * edge = perfectMatching[i];
			// printEdge(vertices, edge);
		}

		return perfectMatching;
	}

	// 4) Combine the edges of M and T to form a connected multigraph H
	//    in which each vertex has even degree.
	// 5) Form an Eulerian circuit in H.
	// 6) Make the circuit found in previous step into a
	//    Hamiltonian circuit by skipping repeated vertices (shortcutting).
	vector<int> createCircuit(vector<City *> vertices, vector<Edge *> mstEdges, vector<Edge *> pmEdges) {
		/////////////////////////////////////////////////////////
		// Based on this algorithm:
		//	http://www.graph-magics.com/articles/euler.php
		// we know graph has 0 odd vertices, so start at any vertex
		// O(V+E) complexity
		/////////////////////////////////////////////////////////

		// cout << "Step 4: Caclulating union of MTS and Perfect Matching" << endl;

		// make copy of original adjlist to use/modify
		unordered_map<int, vector<int> > temp;

		for (int i = 0; i < mstEdges.size(); i++) {
			Edge * edge = mstEdges[i];
			// printEdge(vertices, edge);
			temp[edge->index1()].push_back(edge->index2());
			temp[edge->index2()].push_back(edge->index1());
		}

		for (int i = 0; i < pmEdges.size(); i++) {
			Edge * edge = pmEdges[i];
			// printEdge(vertices, edge);
			temp[edge->index1()].push_back(edge->index2());
			temp[edge->index2()].push_back(edge->index1());
		}

		// cout << "Step 5: Form an Eulerian circuit in H" << endl;

		// Repeat until the current vertex has no more neighbors and the stack is empty.
		int pos = 0;
		vector<int> path;
		stack<int> stk;
		while (!stk.empty() || temp[pos].size() > 0 ) {
			// If current vertex has no neighbors -
			if (temp[pos].size() == 0) {
				// add it to circuit,
				path.push_back(pos);
				// remove the last vertex from the stack and set it as the current one.
				int last = stk.top();
				stk.pop();
				pos = last;
			}
			// Otherwise (in case it has neighbors)
			else {
				// add the vertex to the stack,
				stk.push(pos);
				// take any of its neighbors,
				int neighbor = temp[pos].back();
				// remove the edge between selected neighbor and that vertex,
				temp[pos].pop_back();
		        for (unsigned int i = 0; i < temp[neighbor].size(); i++)
		            if (temp[neighbor][i] == pos) { // find position of neighbor in list
		        	    temp[neighbor].erase (temp[neighbor].begin() + i); // remove it
		                break;
		            }
				// and set that neighbor as the current vertex.
		        pos = neighbor;
			}
		}

		path.push_back(pos);

		// cout << "Step 6: Make Hamiltonian cuircuit" << endl;

		// remove visited nodes from Euler tour
		unordered_map<int, bool> visited;

		int root = path.front();
		vector<int>::iterator curr = path.begin();
		vector<int>::iterator next = path.begin()+1;
		visited[root] = true;

		// loop until the end of the circuit list is reached
		while ( next != path.end() ) {
			// if we haven't been to the next city yet, go there
			if (!visited[*next]) {
				curr = next;
				visited[*curr] = true;
				next = curr + 1;
			} else {
				next = path.erase(next); // remove it
			}
		}

		return path;
	}

public:
	vector<City *> route(vector<City *> cities) {
		if (cities.size() < 3) {
			return cities;
		} 

		vector<Edge *> edges;
		for(int i = 0; i < cities.size(); i++) {
			City * city1 = cities[i];
			for (int j = i; j < cities.size(); j++) {
				if (i != j) {
					City * city2 = cities[j];
					double distance = city1->distance(city2);
					Edge * edge = new Edge(i, j, distance);
					edges.push_back(edge);
				}
			}
		}

		vector<Edge *> mstEdges = kruskalsMst(cities, edges);
		vector<Edge *> pmEdges = perfectMatching(cities, edges, mstEdges);
		vector<int> circuit = createCircuit(cities, mstEdges, pmEdges);

		vector<City *> finalCircuit;
		for (int i = 0; i < circuit.size(); i++) {
			int index = circuit[i];
			City * city = cities[index];
			finalCircuit.push_back(city);
		}

		return finalCircuit;
	}
};

Router * Router::christofidesRouter() {
	return new ChristofidesRouter();
}

