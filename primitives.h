#ifndef PRIMITIVES_H
#define PRIMITIVES_H

#include <math.h>

class Vec2{
public:
	double x,y;

	Vec2(){ x = y = 0.0;	}
	Vec2(double x,double y){ this->x = x; this->y = y; }
	~Vec2(){}

	Vec2 &operator=(const Vec2 &v){ this->x = v.x;this->y = v.y; return *this; }
	Vec2 &operator+=(const Vec2 &v){ this->x += v.x;this->y += v.y; return *this; }
	Vec2 &operator-=(const Vec2 &v){ this->x -= v.x;this->y -= v.y; return *this; }

	Vec2 operator+(const Vec2 &v){ return Vec2(this->x + v.x,this->y + v.y); }
	Vec2 operator-(const Vec2 &v){ return Vec2(this->x - v.x,this->y - v.y); }

	Vec2 times(double a){ return Vec2(this->x * a,this->y * a); }
	double dot(Vec2 &v){ return this->x * v.x + this->y * v.y; }
	double cross(Vec2 &v){ return this->x * v.y - this->y * v.x; }
	double norm1(){ return fabs(x) + fabs(y); }
	double norm2(){ return sqrt(x*x+y*y); }
	double angle(){ return y >= 0 ? acos(x/norm2()) : 2*M_PI - acos(x/norm2()); }
};

class Primitives{
public:
	Vec2 start_point;
	double start_angle;
	double length;

	Primitives(){ length = 0.0; }
	Primitives(Vec2 sp,double sa,double l){
		start_point = sp; start_angle = sa; length = l;
	}
	virtual ~Primitives(){}

	virtual double distance(Vec2 v){};
};

class Line : public Primitives{
public:
	Line() : Primitives(){}
	Line(Vec2 sp,double sa,double l) : Primitives(sp,sa,l){}
	~Line(){}


	double distance(Vec2 v) override{
		Vec2 vv = start_point - v;
		Vec2 ab(length*cos(start_angle),length*sin(start_angle));
		double aabb = ab.x*ab.x + ab.y*ab.y;
		double t = -(ab.x*vv.x + ab.y*vv.y);
		if(t < 0){
			return vv.norm2();	
		}
		else if(t > aabb){
			return  (vv + ab).norm2();
		}
		else{
			return (ab.x*vv.y - ab.y*vv.x)/sqrt(aabb);
		}
	}
};

class Arc : public Primitives{
public:
	double start_curvature; // 正で反時計回り、負で時計回り

	Arc() : Primitives(){ start_curvature = 0.0; }
	Arc(Vec2 sp,double sa,double l,double sc) : Primitives(sp,sa,l) { start_curvature = sc; }
	~Arc(){}
	
	double radius(){ return fabs(1.0/start_curvature); }
	Vec2 centerPos(){
		double r = radius();
		double a = start_curvature > 0 ? start_angle + M_PI/2 : start_angle - M_PI/2;
		return Vec2(start_point.x + r*cos(a),start_point.y + r*sin(a));
	}
	// 始点から円の中心へのベクトル
	Vec2 vecToCenter(){
		double r = radius();
		double a = start_curvature > 0 ? start_angle + M_PI/2 : start_angle - M_PI/2;
		return Vec2(r*cos(a),r*sin(a));
	}

	double distance(Vec2 v) override{
		double r = radius();
		Vec2 vc = vecToCenter();
		Vec2 vcinv = vc.times(-1);
		Vec2 center = start_point + vc;
		double astart = vcinv.angle();
		double arot = length * start_curvature; // 回転角(時計回り)
		if(arot < 0){ // 始点->終点が反時計回りになるよう調整
			astart += arot;
			arot = -arot;
		}

		Vec2 tv = v - center; // 中心が(0,0)に来るように平行移動して考える
		double atv = tv.angle() - astart; // start_pointがx軸上になるよう回転
		atv = atv < 0 ? atv + 2*M_PI : atv; // 0 <= atv <= 2pi
		atv = atv > 2*M_PI ? atv - 2*M_PI : atv;
		if(atv < arot || arot >= 2*M_PI){
			return fabs(tv.norm2()-r);
		}
		else{
			double acenter = arot/2 + M_PI;
			Vec2 vv(tv.norm2() * cos(atv),tv.norm2()* sin(atv));
			if(atv < acenter){ // 終点に近い
				return (vv - Vec2(r*cos(arot),r*sin(arot))).norm2();
			}
			else{ // 始点に近い
				return (vv - Vec2(r,0)).norm2();
			}
		}

	}
};

class Clothoid : public Primitives{
public:
	double start_curvature;
	double end_curvature;

	Clothoid() : Primitives() { start_curvature = end_curvature = 0.0; }
	Clothoid(Vec2 sp,double sa,double l,double sc,double ec) : Primitives(sp,sa,l) { start_curvature = sc; end_curvature = ec; }
	~Clothoid(){}
};

#endif
