#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <queue>
#include <map>
#include <algorithm>
#include <functional>
#include <limits.h>
#include <float.h>
#include "primitives.h"

#define TYPE_NONE 0
#define TYPE_LINE 1
#define TYPE_ARC 2
#define TYPE_CLOTHOID 4

class Vertex{
public:
	Primitives * prim;
	int type;
	double cost;
	int begin;
	int end;
	Vertex(){
		prim = NULL;
		type = TYPE_NONE;
		cost = 0;
	}
	Vertex(Primitives* prim,int type,double cost,int begin,int end){
		this->prim = prim;
		this->type = type;
		this->cost = cost;
		this->begin = begin;
		this->end = end;
	}
};

class Vertice{
public:
	struct PInfo{
		Points pt;
		int begin;
		int end;
		bool contain_corner;
	};
	
	std::vector<Vertex> vs; // 最後二つは始点と終点
	std::vector<PInfo> ps;
	std::vector<Line> lines;
	std::vector<Arc> arcs;
	std::vector<Clothoid> clothoids;

	int size(){
		return vs.size();
	}
	Vertice(){
		vs = std::vector<Vertex>();
		ps = std::vector<PInfo>();
		lines = std::vector<Line>();
		arcs = std::vector<Arc>();
		clothoids = std::vector<Clothoid>();
	}

	Vertice(Points& pt){
		PointsSubset(pt);
		int end = ps.size();
		lines = std::vector<Line>(end);
		arcs = std::vector<Arc>(end);
		clothoids = std::vector<Clothoid>(end);
		vs = std::vector<Vertex>(end*3+2);
//		std::cerr << end << std::endl;
		for(int i = 0; i < end; i++){
			if(ps[i].contain_corner){
				vs[i*3] =Vertex(); 
				vs[i*3+1] =Vertex(); 
				vs[i*3+2] =Vertex(); 
			}
			else{
				double alpha = 1.0;
				lines[i] = fitLine(ps[i].pt);
//				std::cerr << "L" << ps[i].pt.size() << std::endl;
				vs[i*3] = Vertex(&lines[i],TYPE_LINE,alpha*1+lines[i].GetScore(ps[i].pt),ps[i].begin,ps[i].end);
				arcs[i] = fitArc(ps[i].pt);
				vs[i*3+1] = Vertex(&arcs[i],TYPE_ARC,alpha*2+arcs[i].GetScore(ps[i].pt),ps[i].begin,ps[i].end);
//				std::cerr << "SC" << arcs[i].GetScore(ps[i].pt) << std::endl;
//				clothoids[i] = fitClothoid(ps[i].pt);
				clothoids[i] = Clothoid();
				vs[i*3+2] = Vertex(&clothoids[i],TYPE_CLOTHOID,alpha*4+clothoids[i].GetScore(ps[i].pt),ps[i].begin,ps[i].end);
			}
			if(i % 100 == 0)
				std::cerr << i << std::endl;
		}
		vs[vs.size()-2] = Vertex();
		vs[vs.size()-1] = Vertex();
	}
	
	void PointsSubset(Points& pt){
		int end = pt.size();
		ps = std::vector<PInfo>(end*(end-1)/2);
		int cntr = 0;
		for(int i = 0; i < end; i++){
			for(int j = i+2; j <= end; j++){
				pt.trim(i,j,ps[cntr].pt);
				ps[cntr].begin = i;
				ps[cntr].end = j;
				ps[cntr].contain_corner = pt.containCorner(i,j);
				cntr++;
			}
		}
	}
};

typedef std::pair<int,double> Pr; // to,cost
typedef std::vector<Pr> VPr;
typedef std::vector<VPr> VVPr; // adjacent list

class Edges{
public:
	VVPr es; // adjacent list
	
	Edges(){
		es = VVPr();
	}
	Edges(Vertice& vs,Points& pt){
		es = VVPr(vs.vs.size());
		for(int i = 0; i < vs.vs.size()-2; i++){
			es[i] = VPr();
				if(vs.ps[i/3].contain_corner){
					continue; // 途中に角を含む
				}
			for(int j = i+1; j < vs.vs.size()-2; j++){
				if(vs.ps[j/3].contain_corner){
					continue; // 途中に角を含む
				}
				if(vs.vs[i].end - 1 == vs.vs[j].begin && std::find(pt.corner.begin(),pt.corner.end(),vs.vs[j].begin) != pt.corner.end()){
				// 角でつながるとき
				// TODO
					double add = (vs.vs[j].prim->start_point - vs.vs[i].prim->GetEndPos()).norm1();
					es[i].push_back(std::make_pair(j,(vs.vs[i].cost + vs.vs[j].cost)/2+add));
	//				std::cerr << "FE" << vs.vs[j].type << " " << vs.vs[i].type << " " << vs.vs[i].prim->GetEndPos().norm1() << " " << vs.vs[j].prim->start_point.norm1()  << " " << add << " " << (vs.vs[i].cost + vs.vs[j].cost)/2+add << std::endl;
				}
				else if(vs.vs[i].end -1 == vs.vs[j].begin && std::find(pt.corner.begin(),pt.corner.end(),vs.vs[j].begin) == pt.corner.end()){
				// 辺でつながるとき
				// TODO
				//	double add = (fabs(vs.vs[i].prim->length)+fabs(vs.vs[j].prim->length))/2*
				//		((vs.vs[j].prim->start_point - vs.vs[i].prim->GetEndPos()).norm1() +
				//		fabs(vs.vs[j].prim->start_angle - vs.vs[i].prim->GetEndAngle()));
					double add = 
						((vs.vs[j].prim->start_point - vs.vs[i].prim->GetEndPos()).norm1() +
						fabs(vs.vs[j].prim->start_angle - vs.vs[i].prim->GetEndAngle()));
					es[i].push_back(std::make_pair(j,(vs.vs[i].cost + vs.vs[j].cost)/2 + add));
				}
			}
		}

		int idxS = vs.vs.size()-2;
		int idxE = vs.vs.size()-1;
		es[idxS] = VPr();
		es[idxE] = VPr();
		for(int i = 0; i < vs.vs.size()-2; i++){
			if(vs.ps[i/3].contain_corner){
				continue;
			}
			else if(vs.vs[i].begin == 0){
				es[idxS].push_back(std::make_pair(i,1.0)); // 始点と結ぶ
			}
			else if(vs.vs[i].end == pt.size()){ 
				es[i].push_back(std::make_pair(idxE,1.0)); // 終点と結ぶ
			}
		}
	}
};


class Graph{
private:
	typedef std::pair<double,int> P;
public:
	Vertice V;
	Edges E;
	std::vector<int> shortestPath;

	Graph(Points& pt){
		std::cerr << "V" << std::endl;
		V = Vertice(pt);
		std::cerr << "E" << std::endl;
		E = Edges(V,pt);
		shortestPath = std::vector<int>();
	}
	
	void dijkstra(){
		shortestPath = std::vector<int>();
		int s = V.size()-2; // 始点
		int t = V.size()-1; // 終点
		double* d = new double[V.size()];
		int* from = new int[V.size()];
		std::priority_queue<P,std::vector<P>,std::greater<P> > que;
		for(int i = 0; i < V.size(); i++){
			d[i] = DBL_MAX;
			from[i] = INT_MAX;
		}
		d[s] = 0;
		que.push(P(0,s));

		std::cerr << "dijkstra" << std::endl;

		while(!que.empty()){
			P p = que.top(); que.pop();
			int v = p.second;
			if(d[v] < p.first) continue;
			for(int i = 0; i < E.es[v].size(); i++){
				Pr e = E.es[v][i]; // <int,double>
	//			std::cerr << v << " " << d[e.first] <<  " " << d[v] << " " << e.second << " " << E.es[v].size() << std::endl;
				if(d[e.first] > d[v] + e.second){
					d[e.first] = d[v] + e.second;
					from[e.first] = v;
					que.push(P(d[e.first],e.first));
				}
			}
		}

//		for(int i = 0; i < V.size(); i++){
//			std::cerr <<i << " " << d[i] << " " << from[i] << " ";
//			std::cerr << std::endl;
//		}
		std::cerr << "back: s,t =" << s << "," << t << std::endl;
		
		int back = t;
		while(back != s){
			shortestPath.push_back(back);
			if(back > t){
				std::runtime_error("dijkstra err");
			}
			std::cerr << back << std::endl;
			back = from[back];
		}
		shortestPath.push_back(s);
		std::reverse(shortestPath.begin(),shortestPath.end());
		delete[] d;
		delete[] from;
	}

	Points GetAllPoints(){
		dijkstra();
		int n = 50;
		Points res((shortestPath.size()-2)*n+1);
		for(int i = 1; i < shortestPath.size()-1; i++){ // 始点と終点を除く	
			Points pt;
			V.vs[shortestPath[i]].prim->GetPoints(n,pt);
			for(int j = 0; j < n; j++){
				res.at((i-1)*n+j) = pt.at(j);
			}
		}
		res.at(res.size()-1) = res.at(0);
		return res;
	}
};


#endif
