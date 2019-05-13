#include <bits/stdc++.h>
#include "primitives.h"
#include "graph.h"
using namespace std;

int main(){
//	FILE* fp1 = fopen("./debug/boundary.txt","r");
//	FILE* fp2 = fopen("./debug/corners.txt","r");
//	FILE* fp1 = fopen("./16/shape3_boundary.txt","r");
//	FILE* fp2 = fopen("./16/shape3_corners.txt","r");
	FILE* fp1 = fopen("./fighter/32_boundary.txt","r");
	FILE* fp2 = fopen("./fighter/32_corners.txt","r");
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
	cerr << "Points " << endl;
	pl.print();
	Graph G(pl);
	cerr << "C " << pt.corner.size() << endl;
	cerr << "Path" << endl;
	Points res = G.GetAllPoints();
	if(res.size() >= pl.size()){
		for(int i = 0; i < res.size(); i++){
			cout << res.at(i).x << " " << res.at(i).y << " ";
			if(i < pl.size())
				cout << pl.at(i).x << " " << pl.at(i).y << endl;
			else
				cout << endl;
		}
	}
	else{
		for(int i = 0; i < pl.size(); i++){
			cout << pl.at(i).x << " " << pl.at(i).y << " ";
			if(i < res.size())
				cout << res.at(i).x << " " << res.at(i).y << endl;
			else
				cout << endl;
		}
	}

	int n = 40;
	/*
	Points giv(5);
	giv.at(0) = Vec2(0,0);
	giv.at(1) = Vec2(1,0);
	giv.at(2) = Vec2(2,1);
	giv.at(3) = Vec2(2,2);
	giv.at(4) = Vec2(2,3);
	Clothoid c(Vec2(3.0,2.0),0.78,1,-1.0,1.0);
	Points pt;
	c.GetPoints(n,pt);
	Clothoid ln = fitClothoid(pt);
	Points pt2;
	ln.GetPoints(n,pt2);
	for(int i = 0; i < n; i++){
		cout << pt2.at(i).x << " " << pt2.at(i).y << " ";
		if(i < pt.size())
			cout << pt.at(i).x << " " << pt.at(i).y << endl;
		else
			cout << endl;
//		cout << pt.at(i).x << " " << pt.at(i).y << " " << endl;
	}
//	cout << c.distance(Vec2(9.5,2.5)) << endl;
*/
	return 0;
}
