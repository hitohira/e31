#include <bits/stdc++.h>
#include "primitives.h"
using namespace std;

int main(){
	int n = 100;
	Clothoid c(Vec2(3.0,2.0),0.78,1,-0.2,2.0);
	Points pt = c.GetPoints(n);
	for(int i = 0; i < n; i++){
		cout << pt.at(i).x << " " << pt.at(i).y << endl;
	}
	return 0;
}
