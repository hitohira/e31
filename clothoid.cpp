#include <bits/stdc++.h>
#include "primitives.h"
using namespace std;

int main(){
	Clothoid c(Vec2(3.0,2.0),1.4,2.0,0.0,4.0);
	for(int i = -100; i < 100; i++){
		Vec2 v = c.posByT(i/20.0);
		cout << v.x << " " << v.y << endl;
	}
	return 0;
}
