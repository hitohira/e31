#include <bits/stdc++.h>
#include "primitives.h"
using namespace std;

int main(){
	Vec2 v1;
	Vec2 v2(2.0,4.0);
 	Vec2 v = (v1 += v2);
	v2 = v;
	cout << v2.norm2() << endl;
	Primitives *p = new Arc(v2,M_PI/2,M_PI/2,-1);
	cout << p->length << " " << p->start_point.norm2() << endl;
	double d = p->distance(Vec2(2,5));
	cout << d << endl;
	delete p;
	return 0;
}
