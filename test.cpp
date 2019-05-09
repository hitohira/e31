#include <bits/stdc++.h>
#include "primitives.h"
using namespace std;

int main(){
	int n = 40;
	Clothoid c(Vec2(3.0,2.0),0.78,1,-8.0,1.0);
	Points pt;
	c.GetPoints(n,pt);
	Clothoid ln = fitClothoid(pt);
	Points pt2;
	ln.GetPoints(n,pt2);
	for(int i = 0; i < n; i++){
		cout << pt.at(i).x << " " << pt.at(i).y << " ";
		cout << pt2.at(i).x << " " << pt2.at(i).y << endl;
	}
//	cout << c.distance(Vec2(9.5,2.5)) << endl;
	return 0;
}
