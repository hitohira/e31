#ifndef PRIMITIVES_H
#define PRIMITIVES_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include <stdexcept>
#include <vector>
#include "Eigen/Core"
#include "Eigen/LU"

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
	Vec2 rot(double theta){ return Vec2(x*cos(theta) - y*sin(theta),x*sin(theta)+y*cos(theta)); }
};


class Points{
private:
	int sz;
	Vec2 * idx;

public:
	std::vector<int> corner;

	Points(){ sz = 0; idx = NULL; corner = std::vector<int>(); }
	Points(int n){ sz = n; idx = (Vec2*)malloc(sizeof(Vec2)*n); }
	~Points(){ 
		if(idx != NULL){ 
			free(idx);
			idx = NULL; 
		} 
	}
	
	Vec2& at(unsigned int i){
		if(i >= sz){
			throw std::runtime_error("invalid access");
		}
		return idx[i];
	}
	Vec2 at2(unsigned int i) const {
		if(i >= sz){
			throw std::runtime_error("invalid access");
		}
		return idx[i];
	}
	int size(){
		return sz;
	}

	Points &operator=(const Points &p)
	{ 
		if(idx != NULL){
			free(idx);
			idx = NULL;
		}
		this->sz = p.sz;
		this->idx = (Vec2*)malloc(sizeof(Vec2)*p.sz);
		for(int i = 0; i < this->sz; i++){
			at(i) = p.at2(i);
		}
		return *this;
	}

	Vec2 average(){
		Vec2 sum;
		for(int i = 0; i < sz; i++){
			sum += idx[i];
		}
		return sum.times(1.0/sz);
	}
	
	int readPoints(FILE* fp){
		double d1,d2;
		std::vector<double> v1 = std::vector<double>();
		std::vector<double> v2 = std::vector<double>();
		while(fscanf(fp,"%lf%lf",&d1,&d2) == 2){
			v1.push_back(d1);
			v2.push_back(d2);
		}
		if(idx != NULL){
			free(idx);
		}
		sz = v1.size();
		idx = (Vec2*)malloc(sizeof(Vec2)*sz);
		for(int i = 0; i < sz; i++){
			idx[i] = Vec2(v1[i],v2[i]);
		}
		return sz;
	}
	int readCornerIdx(FILE* fp){
		int p;
		if(!corner.empty()){
			std::vector<int>().swap(corner);
		}
		while(fscanf(fp,"%d",&p) == 1){
			corner.push_back(p);
		}
		return corner.size();
	}
	
	// 先頭が角になるように並びかえ、先頭と最後が同じ点になるようにする
	Points rearrenge(){
		Points pt(sz+1);
		if(corner.empty()){
			for(int i = 0; i < sz; i++){
				pt.at(i) = at(i);
				pt.at(sz) = at(0);
			}
		}
		else{
			int top = corner[0];
			for(int i = top; i < sz; i++){
				pt.at(i-top) = at(i);
			}
			for(int i = 0; i < top; i++){
				pt.at(i+sz-top) = at(i);
			}
			pt.at(sz) = at(top);
			for(int i = 0; i < corner.size(); i++){
				pt.corner.push_back(corner[i]-top);
			}
			pt.corner.push_back(sz);
		}
		return pt;
	}
	// 部分切り出し
	void trim(int begin,int end,Points &dst){
		dst = Points(end-begin);
		for(int i = begin; i < end; i++){
			dst.at(i-begin) = at(i);
		}
	}

	bool containCorner(int begin,int end){
		for(int i = 0; i < corner.size(); i++){
			if(corner[i] > begin && corner[i] < end-1){
				return true;
			}
		}
		return false;
	}

};


class Primitives {
public:
	Vec2 start_point;
	double start_angle;
	double length;

	Primitives(){ start_point = Vec2(0,0);start_angle = 0.0;length = 0.0; }
	Primitives(Vec2 sp,double sa,double l){
		start_point = sp; start_angle = sa; length = l;
	}
	virtual ~Primitives(){}

	virtual double distance(Vec2 v){ return 0.0; }

	virtual void GetPoints(int n,Points& pt){ pt = Points();}

	virtual double GetScore(Points& pt) { return 0.0; }

	virtual double GetEndAngle(){ return 0.0 ;}

	virtual Vec2 GetEndPos(){ return Vec2(); }
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
			return fabs((ab.x*vv.y - ab.y*vv.x))/sqrt(aabb);
		}
	}
	
	void GetPoints(int n,Points& pt) override{
		pt = Points(n);
		pt.at(0) = start_point;
		for(int i = 1; i < n; i++){
			double l = length / (n-1) * i;
			pt.at(i) = start_point + Vec2(l*cos(start_angle),l*sin(start_angle));
		}
	}

	double GetScore(Points& pt) override{
		double sm = 0.0;
		for(int i = 0; i < pt.size(); i++){
			sm += distance(pt.at(i));
		}
		return sm;
	}
	double GetEndAngle() override{
		return start_angle;
	}
	Vec2 GetEndPos() override{
		return start_point + Vec2(length*cos(start_angle),length*sin(start_angle));
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
		if(fabs(radius()) < 1e-10){
			Line ln(start_point,start_angle,length);
			return ln.distance(v);
		}
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

	void GetPoints(int n, Points& pt) override{
		if(fabs(radius()) < 1e-10){
			Line ln(start_point,start_angle,length);
			return ln.GetPoints(n,pt);
		}
		pt = Points(n);
		pt.at(0) = start_point;
		double r = radius();
		Vec2 vc = vecToCenter();
		Vec2 vcinv = vc.times(-1);
		Vec2 center = start_point + vc;
		double arot = length * start_curvature; // 回転角(時計回り)
		for(int i = 1; i < n; i++){
			double a = arot / (n-1) * i;
			pt.at(i) = center + vcinv.rot(a);	
		}
	}
	double GetScore(Points& pt) override {
		double sm = 0.0;
		for(int i = 0; i < pt.size(); i++){
			sm += distance(pt.at(i));
		}
		return sm;
	}

	double GetEndAngle() override{
		double theta = length / radius();
		double angle;
		if(start_curvature >= 0.0){
			angle = start_angle + theta;
		}
		else{
			angle = start_angle - theta;
		}
		if(angle > 2*M_PI){
			angle -= 2*M_PI;
		}
		else if(angle < 0.0){
			angle += 2*M_PI;
		}
		return angle;
	}
	Vec2 GetEndPos() override{
		if(fabs(radius()) < 1e-10){
			Line ln(start_point,start_angle,length);
			return ln.GetEndPos();
		}
		double r = radius();
		Vec2 vc = vecToCenter();
		Vec2 vcinv = vc.times(-1);
		Vec2 center = start_point + vc;
		double arot = length * start_curvature; // 回転角(時計回り)
		return center + vcinv.rot(arot);	
	}
};

class Clothoid : public Primitives{
public:
	double start_curvature;
	double end_curvature;

	Clothoid() : Primitives() { start_curvature = end_curvature = 0.0; }
	Clothoid(Vec2 sp,double sa,double l,double sc,double ec) : Primitives(sp,sa,l) { start_curvature = sc; end_curvature = ec; }
	~Clothoid(){}

	Vec2 posByT(double t){
		double B2 = fabs(length / (M_PI*(end_curvature-start_curvature)));
		double B = sqrt(B2);
		double start_t = start_curvature * B;
		double end_t =  end_curvature * B;

		double len = start_curvature <= end_curvature ? length : -length;
		double start_R = 1.0 / start_curvature;
		double end_R = 1.0 / end_curvature;
		double start_L = end_R * len / (start_R - end_R);
		double start_tau = start_L / (2*start_R);
		double rota = start_tau - start_angle; 

		Vec2 v = normalPosByBT(B,t) - normalPosByBT(B,start_t);
		v = v.rot(-rota);
		v = v + start_point;
		return v;
	}
	
	Vec2 normalPosByT(double t){
		double B2 = fabs(length / (M_PI*(end_curvature-start_curvature)));
		double B = sqrt(B2);
		return normalPosByBT(B,t);
	}

	Vec2 normalPosByBT(double B,double t){
		bool sgn = t >= 0;
		t = fabs(t);
		double t2 = t*t;
		double t3 = t2*t;
		double R = (0.506*t+1)/(1.79*t2+2.054*t+sqrt(2.0));
		double A = 1.0 / (0.803*t3+1.886*t2+2.524*t+2);
		double C = 0.5 - R*sin(0.5*M_PI*(A-t2));
		double S = 0.5 - R*cos(0.5*M_PI*(A-t2));
		double x = M_PI*B*C;
		double y = M_PI*B*S;
		return  sgn ? Vec2(x,y) : Vec2(-x,-y);
	}
	Vec2 derByBT(double B,double t){
		double mtt2 = M_PI/2*t*t;
		double mb = M_PI*B;
		return Vec2(mb*cos(mtt2),mb*sin(mtt2));
	}
	Vec2 der2ByBT(double B,double t){
		double ptt2 = M_PI/2*t*t;
		double ppbt = M_PI*M_PI*B*t;
		return Vec2(-ppbt*sin(ptt2),ppbt*cos(ptt2));
	}

	double f(Vec2 p,double B,double t){
		Vec2 v = normalPosByBT(B,t);
		Vec2 tmp = v-p;
		return tmp.x*tmp.x + tmp.y*tmp.y;
	}
	double f1(Vec2 p, double B, double t){
		Vec2 v = normalPosByBT(B,t);
		Vec2 d = derByBT(B,t);
		return -2*(p.x-v.x)*d.x -2*(p.y-v.y)*d.y;
	}
	double f2(Vec2 p,double B,double t){
		Vec2 v = normalPosByBT(B,t);
		Vec2 d2 = der2ByBT(B,t);
		return -2*(p.x-v.x)*d2.x - 2*(p.y-v.y)*d2.y + 2*M_PI*M_PI*B*B;
	}

	double dist_newton(Vec2 p,double B,double min_t,double max_t,double t){
		double val = f1(p,B,t);
		int cntr = 0;
		while(fabs(val) > 1e-4 && cntr < 20){
			t -= val / f2(p,B,t);
			val = f1(p,B,t);
			cntr++;
		}
		if(t < min_t){
			return sqrt(f(p,B,min_t));
		}
		if(t > max_t){
			return sqrt(f(p,B,max_t));
		}
		return sqrt(f(p,B,t));
	}

	double distance(Vec2 v) override{
		if(fabs(end_curvature - start_curvature) < 1e-10){
			Arc arc(start_point,start_angle,length,start_curvature);
			return arc.distance(v);
		}
		//　パラメータ計算
		double B2 = fabs(length / (M_PI*(end_curvature-start_curvature)));
		double B = sqrt(B2);
		double start_t = start_curvature * B;
		double end_t =  end_curvature * B;

		double len = start_curvature <= end_curvature ? length : -length;
		double start_R = 1.0 / start_curvature;
		double end_R = 1.0 / end_curvature;
		double start_L = end_R * len / (start_R - end_R);
		double start_tau = start_L / (2*start_R);
		double rota = start_tau - start_angle; // 正規形にするために回転させる角度
		// 正規形へ座標変換
		v = v - start_point;
		v = v.rot(rota);
		v = v + normalPosByBT(B,start_t);
		// 距離をiterativeに計算
		double min_t = fmin(start_t,end_t);
		double max_t = fmax(start_t,end_t);
		double d1 = dist_newton(v,B,min_t,max_t,min_t);
		double d2 = dist_newton(v,B,min_t,max_t,max_t);
		double d3 = dist_newton(v,B,min_t,max_t,(min_t+max_t)/2);
		return fmin(d1,fmin(d2,d3));
	}

	void GetPoints(int n,Points& pt)override{
		if(fabs(end_curvature - start_curvature) < 1e-10){
			Arc arc(start_point,start_angle,length,start_curvature);
			return arc.GetPoints(n,pt);
		}
		pt = Points(n);
		pt.at(0) = start_point;
		double B2 = fabs(length / (M_PI*(end_curvature-start_curvature)));
		double B = sqrt(B2);
		double start_t = start_curvature * B;
		double end_t =  end_curvature * B;
		for(int i = 1; i < n; i++){
			double t = start_t + (end_t - start_t) / (n-1) * i;
			pt.at(i) = posByT(t); 
		}
	}
	double GetScore(Points& pt)override{
		double sm = 0.0;
		for(int i = 0; i < pt.size(); i++){
			sm += distance(pt.at(i));
		}
		return sm;
	}
	double GetEndAngle()override{
		double len = start_curvature <= end_curvature ? length : -length;
		double start_R = 1.0 / start_curvature;
		double end_R = 1.0 / end_curvature;
		double start_L = end_R * len / (start_R - end_R);
		double end_L = start_L + len;
		double start_tau = start_L / (2*start_R);
		double rota = start_tau - start_angle; // 正規形にするために回転させる角度
		double angle = start_L / (2*start_R) - rota;
		if(angle < 2*M_PI){
			angle += 2*M_PI;
		}
		else if(angle < 0){
			angle += 2*M_PI;
		}
		return angle;
	}
	Vec2 GetEndPos()override{
		if(fabs(end_curvature - start_curvature) < 1e-10){
			Arc arc(start_point,start_angle,length,start_curvature);
			return arc.GetEndPos();
		}
		double B2 = fabs(length / (M_PI*(end_curvature-start_curvature)));
		double B = sqrt(B2);
		double end_t =  end_curvature * B;
		return posByT(end_t); 
	}
};

Vec2 solveQuadEq(double a,double b,double c){
	double d = b*b-4*a*c;
	if(b > 0){
		double x =(-b-sqrt(d))/(2*a);
		return Vec2(x,c/a/x);
	}
	else{
		double x = (-b+sqrt(d))/(2*a);
		return Vec2(x,c/a/x);
	}
}

Line fitLine(Points &pt){
	Vec2 mbar = pt.average();
	double a = 0;
	double b = 0;
	double c = 0;
	double d = 0;
	// 共分散行列
	for(int i = 0; i < pt.size(); i++){
		Vec2 sub = pt.at(i) -mbar;
		a += sub.x * sub.x;
		b += sub.x * sub.y;
		c += sub.x * sub.y;
		d += sub.y * sub.y;
	}
	Vec2 x = solveQuadEq(1.0,-(a+b),a*d-b*c);
	double l = fmax(x.x,x.y); // 最大固有値
	Vec2 lv(b,-a+l);
	lv = lv.times(1.0/lv.norm2()); // これが傾き
	// mbarをとおる
	Vec2 m1 = pt.at(0) - mbar;
	Vec2 startPos = mbar + lv.times(m1.dot(lv));
	Vec2 mn = pt.at(pt.size()-1) - mbar;
	Vec2 endPos = mbar + lv.times(mn.dot(lv));
	double len = (startPos - endPos).norm2();
	double angle = lv.angle();
	return Line(startPos,angle,len);
}

void derivArc(Arc& arc,Points& pt,Eigen::VectorXd& F,Eigen::MatrixXd& dF,Eigen::MatrixXd& tdF){
	int n = pt.size();
	F = Eigen::VectorXd::Zero(n);
	dF = Eigen::MatrixXd::Zero(n,5);
	tdF = Eigen::MatrixXd::Zero(5,n);

	double dx = 0.00001;
	double dy = 0.00001;
	double da = 0.00001;
	double dl = 0.00001;
	double dc = 0.00001;
	Arc arc_x(Vec2(arc.start_point.x+dx,arc.start_point.y),arc.start_angle,arc.length,arc.start_curvature);
	Arc arc_y(Vec2(arc.start_point.x,arc.start_point.y+dy),arc.start_angle,arc.length,arc.start_curvature);
	Arc arc_a(Vec2(arc.start_point.x,arc.start_point.y),arc.start_angle+da,arc.length,arc.start_curvature);
	Arc arc_l(Vec2(arc.start_point.x,arc.start_point.y),arc.start_angle,arc.length+dl,arc.start_curvature);
	Arc arc_c(Vec2(arc.start_point.x,arc.start_point.y),arc.start_angle,arc.length,arc.start_curvature+dc);
	for(int i = 0; i < n; i++){
		double dist = arc.distance(pt.at(i));
		F(i) = dist;
		// 差分法
		// point.x
		double dist_x = arc_x.distance(pt.at(i));
		dF(i,0) = (dist_x - dist) / dx;
		tdF(0,i) = (dist_x - dist) / dx;
		// point.y
		double dist_y = arc_y.distance(pt.at(i));
		dF(i,1) = (dist_y - dist) / dy;
		tdF(1,i) = (dist_y - dist) / dy;
		// angle
		double dist_a = arc_a.distance(pt.at(i));
		dF(i,2) = (dist_a - dist) / da;
		tdF(2,i) = (dist_a - dist) / da;
		// length
		double dist_l = arc_l.distance(pt.at(i));
		dF(i,3) = (dist_l - dist) / dl;
		tdF(3,i) = (dist_l - dist) / dl;
		// curvature
		double dist_c = arc_c.distance(pt.at(i));
		dF(i,4) = (dist_c - dist) / dc;
		tdF(4,i) = (dist_c - dist) / dc;
	}
}

Vec2 arcCenterBy3points(Vec2 v1,Vec2 v2,Vec2 v3){
	Eigen::Matrix2d A;
	A << 2.0 * (v1.x-v2.x) , 2.0 * (v1.y-v2.y) ,
	     2.0 * (v1.x-v3.x) , 2.0 * (v1.y-v3.y);
	Eigen::Vector2d b;
	b << (v1.x*v1.x + v1.y*v1.y) - (v2.x*v2.x + v2.y*v2.y), 
	     (v1.x*v1.x + v1.y*v1.y) - (v3.x*v3.x + v3.y*v3.y);
	Eigen::PartialPivLU<Eigen::Matrix2d> lu(A);
	Eigen::Vector2d x = lu.solve(b);
	return Vec2(x(0),x(1));	    
}

Arc initArc(Points & pt){
	int n = pt.size();
	if(n < 3){
		Vec2 dir = pt.at(n-1) - pt.at(0);
		Vec2 sp = pt.at(0);
		double r = dir.norm2();
		double len = r * (M_PI / 3.0);
		double angle = dir.rot(-M_PI/6.0).angle();
		return Arc(sp,angle,len,1.0/r);
	}
	else{
		Vec2 v1 = pt.at(0);
		Vec2 v2 = pt.at(n/2);
		Vec2 v3 = pt.at(n-1);
		Vec2 center = arcCenterBy3points(v1,v2,v3);
		Vec2 vtc = center - v1;
		Vec2 v31 = v3-v1;
		Vec2 v21 = v2-v1;
		double sgn = v21.cross(v31) >= 0 ? 1.0 : -1.0; // 曲率の符号
		double r = vtc.norm2();
		double angle = sgn >= 0 ? vtc.rot(-M_PI/2).angle() : vtc.rot(M_PI/2).angle();
		Vec2 stv = v1 - center;
		Vec2 edv = v3 - center;
		double cost = stv.dot(edv) / (stv.norm2() * edv.norm2());
		double sint = sgn * stv.cross(edv) / (stv.norm2() * edv.norm2());
		double theta = sint >= 0 ? acos(cost) : 2*M_PI - acos(cost); // 中心角
		double len = r * theta;
		return Arc(v1,angle,len,sgn/r);
	}
}

Arc fitArc(Points &pt){
	Arc arc = initArc(pt); // 適切な初期値
	Eigen::VectorXd F;
	Eigen::MatrixXd dF; // 要素数,パラメタ数
	Eigen::MatrixXd tdF; // 要素数,パラメタ数
	Eigen::MatrixXd lambdaI = Eigen::MatrixXd::Identity(5,5) * 0.000001; // 時間あればMarquardt-Levenberg method
	double score = arc.GetScore(pt);
//	std::cerr << score << std::endl;
	for(int i = 0; i < 10; i++){
		derivArc(arc,pt,F,dF,tdF);
		Eigen::MatrixXd A = tdF * dF + lambdaI;
		Eigen::VectorXd b = -tdF * F;
		Eigen::PartialPivLU<Eigen::MatrixXd> lu(A);
		Eigen::VectorXd dx = lu.solve(b);
		if(dx.norm() < 1e-4){
			break;
		}
		Arc next(arc.start_point + Vec2(dx(0),dx(1)),arc.start_angle+dx(2),arc.length+dx(3),arc.start_curvature+dx(4));
		int cntr = 0;
		while(score < next.GetScore(pt) && cntr < 10){
			dx *= 0.8;
			next = Arc(arc.start_point + Vec2(dx(0),dx(1)),arc.start_angle+dx(2),arc.length+dx(3),arc.start_curvature+dx(4));
			cntr++;
		}
		if(cntr == 10){
			break;
		}
		score = next.GetScore(pt);
		arc = next;
//		std::cerr << score << std::endl;
	}
	// 最後に適切にクランプする必要あり
	double r = fabs(1.0 / arc.start_curvature);
	Vec2 center = arc.centerPos();
	Vec2 stv = pt.at(0) - center;
	Vec2 edv = pt.at(pt.size()-1) - center;
	Vec2 startPos = center + stv.times(r/stv.norm2());
	Vec2 endPos = center + edv.times(r/edv.norm2());
	Vec2 vtc = center - startPos;
	double angle = arc.start_curvature >= 0.0 ? vtc.rot(-M_PI/2).angle() : vtc.rot(M_PI/2).angle();
	double sgn = arc.start_curvature >= 0.0 ? 1.0 : -1.0;
	double cost = stv.dot(edv) / (stv.norm2() * edv.norm2());
	double sint = sgn * stv.cross(edv) / (stv.norm2() * edv.norm2());
	double theta = sint >= 0 ? acos(cost) : 2*M_PI - acos(cost); // 中心角
	double len = r * theta;
	return Arc(startPos,angle,len,arc.start_curvature);
}

Clothoid initClothoid(Points& pt){
	Arc arc = initArc(pt);
	return Clothoid(arc.start_point-Vec2(0.05*cos(arc.start_angle),0.05*cos(arc.start_angle)),arc.start_angle,arc.length+1,arc.start_curvature,arc.start_curvature*0.999999);
}

void derivClothoid(Clothoid& cl,Points& pt,Eigen::VectorXd& F,Eigen::MatrixXd& dF,Eigen::MatrixXd& tdF){
	int n = pt.size();
	F = Eigen::VectorXd::Zero(n);
	dF = Eigen::MatrixXd::Zero(n,6);
	tdF = Eigen::MatrixXd::Zero(6,n);

	double dx = 0.0001;
	double dy = 0.0001;
	double da = 0.0001;
	double dl = 0.0001;
	double dsc = 0.0001;
	double dec = 0.0001;
	Clothoid c_x(Vec2(cl.start_point.x+dx,cl.start_point.y),cl.start_angle,cl.length,cl.start_curvature,cl.end_curvature);
	Clothoid c_y(Vec2(cl.start_point.x,cl.start_point.y+dy),cl.start_angle,cl.length,cl.start_curvature,cl.end_curvature);
	Clothoid c_a(cl.start_point,cl.start_angle+da,cl.length,cl.start_curvature,cl.end_curvature);
	Clothoid c_l(cl.start_point,cl.start_angle,cl.length+dl,cl.start_curvature,cl.end_curvature);
	Clothoid c_sc(cl.start_point,cl.start_angle,cl.length,cl.start_curvature+dsc,cl.end_curvature);
	Clothoid c_ec(cl.start_point,cl.start_angle,cl.length,cl.start_curvature,cl.end_curvature+dec);
	for(int i = 0; i < n; i++){
		double dist = cl.distance(pt.at(i));
		F(i) = dist;
		// 差分法
		// point.x
		double dist_x = c_x.distance(pt.at(i));
		dF(i,0) = (dist_x - dist) / dx;
		tdF(0,i) = (dist_x - dist) / dx;
		// point.y
		double dist_y = c_y.distance(pt.at(i));
		dF(i,1) = (dist_y - dist) / dy;
		tdF(1,i) = (dist_y - dist) / dy;
		// angle
		double dist_a = c_a.distance(pt.at(i));
		dF(i,2) = (dist_a - dist) / da;
		tdF(2,i) = (dist_a - dist) / da;
		// length
		double dist_l = c_l.distance(pt.at(i));
		dF(i,3) = (dist_l - dist) / dl;
		tdF(3,i) = (dist_l - dist) / dl;
		// curvature
		double dist_sc = c_sc.distance(pt.at(i));
		dF(i,4) = (dist_sc - dist) / dsc;
		tdF(4,i) = (dist_sc - dist) / dsc;
		// curvature
		double dist_ec = c_ec.distance(pt.at(i));
		dF(i,5) = (dist_ec - dist) / dec;
		tdF(5,i) = (dist_ec - dist) / dec;
	}
}

Clothoid fitClothoid(Points& pt){
	Clothoid cl = initClothoid(pt);
	Eigen::VectorXd F;
	Eigen::MatrixXd dF; // 要素数,パラメタ数
	Eigen::MatrixXd tdF; // 要素数,パラメタ数
	Eigen::MatrixXd lambdaI = Eigen::MatrixXd::Identity(6,6) * 0.000001; // 時間あればMarquardt-Levenberg method
	double score = cl.GetScore(pt);
	for(int i = 0; i < 10; i++){
		derivClothoid(cl,pt,F,dF,tdF);
		Eigen::MatrixXd A = tdF * dF + lambdaI;
		Eigen::VectorXd b = -tdF * F;
		Eigen::PartialPivLU<Eigen::MatrixXd> lu(A);
		Eigen::VectorXd dx = lu.solve(b);
		if(dx.norm() < 1e-4){
			break;
		}
		Clothoid next = Clothoid(cl.start_point + Vec2(dx(0),dx(1)),cl.start_angle+dx(2),
		              cl.length+dx(3),cl.start_curvature+dx(4),cl.end_curvature+dx(5));
		int cntr = 0;
		while(score < next.GetScore(pt) && cntr < 10){
			dx *= 0.8;
			next = Clothoid(cl.start_point + Vec2(dx(0),dx(1)),cl.start_angle+dx(2),
		  		cl.length+dx(3),cl.start_curvature+dx(4),cl.end_curvature+dx(5));
		}
		if(cntr == 10){
			break;
		}
		score = next.GetScore(pt);
		cl = next;
//		std::cerr << score << std::endl;
	}
	// TODO クランプ
	

	return cl;
}


#endif
