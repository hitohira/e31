#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include "primitives.h"

#define TYPE_LINE 1
#define TYPE_ARC 2
#define TYPE_CLOTHOID 4

class Vertex{
public:
	Primitives * prim;
	int type;
	double cost;
	int start;
	int end;
	Vertex(Primitives* prim,int type,double cost,int start,int end){
		this->prim = prim;
		this->type = type;
		this->cost = cost;
		this->start = start;
		this->end = end;
	}
};

class Vertice{
public:
	std::vector<Vertex> vs;

	Vertice(){
		vw = std::vector<Vertex>();
	}
};

class Edges{
public:
	std::vector<std::vector<double> > es; // adjacent list

	Edges(){
		es = std::vector<std::vector<double> >();
	}
};

class Graph{
public:
	Vertice V;
	Edges E;

	Graph(Vertice V){
		this->V = V;
		E = Edges();
	}
};


#endif
