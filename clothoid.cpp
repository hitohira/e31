#include <bits/stdc++.h>
#include "primitives.h"
using namespace std;

int main(){
	int n = 40;
	Clothoid c(Vec2(3.0,2.0),0.78,1,-0.2,2.0);
	Points pt = c.GetPoints(n);
	Line ln = fitLine(pt);
	Points pt2 = ln.GetPoints(n);
	for(int i = 0; i < n; i++){
		cout << pt.at(i).x << " " << pt.at(i).y << " ";
		cout << pt2.at(i).x << " " << pt2.at(i).y << endl;
	}
	return 0;
}
