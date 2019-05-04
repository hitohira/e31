#include <bits/stdc++.h>
#include "primitives.h"
using namespace std;

int main(){
	Vec2 v1;
	Vec2 v2(2.0,4.0);
	Vec2 v3 = v2.rot(1.3);
 	Vec2 v = (v1 += v2);
	v2 = v;
	cout << v2.norm2() << " " << v2.x << " " << v3.x << endl;
	Primitives *p = new Line(v2,M_PI/4,sqrt(2));
	cout << p->length << " " << p->start_point.norm2() << endl;
	double d = p->distance(Vec2(2,5));
	cout << d << endl;
	delete p;

	Points pt(4);
	pt.at(1) = v1;
	pt.at(3) = v2;
	cout << pt.at(1).x << endl;
	return 0;
}
