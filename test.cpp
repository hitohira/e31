#include <bits/stdc++.h>
#include "primitives.h"
#include "graph.h"
using namespace std;

int main(){
	FILE* fp1 = fopen("./debug/boundary.txt","r");
	FILE* fp2 = fopen("./debug/corners.txt","r");
	if(fp1 == NULL || fp2 == NULL){
		cerr << "fopen failed" << endl;
		return -1;
	}
	Points pt;
	pt.readPoints(fp1);
	pt.readCornerIdx(fp2);
	fclose(fp1);
	fclose(fp2);


	cerr << "G" << endl;
	Points pl = pt.rearrenge();
	Graph G(pl);
	cerr << "D" << endl;
	G.dijkstra();
	cerr << "Path" << endl;
	for(int i = 0; i < G.shortestPath.size(); i++){
		cerr << G.shortestPath[i] << endl;
	}
	/*
	int n = 40;
	Points giv(5);
	giv.at(0) = Vec2(0,0);
	giv.at(1) = Vec2(1,0);
	giv.at(2) = Vec2(2,1);
	giv.at(3) = Vec2(2,2);
	giv.at(4) = Vec2(2,3);
	Clothoid c(Vec2(3.0,2.0),0.78,1,-1.0,1.0);
	Points pt;
	c.GetPoints(n,pt);
	Arc ln = fitArc(giv);
	Points pt2;
	ln.GetPoints(n,pt2);
	for(int i = 0; i < n; i++){
		cout << pt2.at(i).x << " " << pt2.at(i).y << " ";
		if(i < 5)
			cout << giv.at(i).x << " " << giv.at(i).y << endl;
		else
			cout << endl;
//		cout << pt.at(i).x << " " << pt.at(i).y << " " << endl;
	}
//	cout << c.distance(Vec2(9.5,2.5)) << endl;
*/
	return 0;
}
