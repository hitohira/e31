#include <bits/stdc++.h>
#include "primitives.h"
#include "graph.h"
using namespace std;

int main(int argc,char** argv){
	// ファイル取得
//	FILE* fp1 = fopen("./debug/boundary.txt","r");
//	FILE* fp2 = fopen("./debug/corners.txt","r");
//	FILE* fp1 = fopen("./16/shape3_boundary.txt","r");
//	FILE* fp2 = fopen("./16/shape3_corners.txt","r");
	FILE* fp1;
	FILE* fp2;
	if(argc == 3){
		cerr << "using ext file" << endl;
		cerr << argv[1] << endl;
		cerr << argv[2] << endl;
		fp1 = fopen(argv[1],"r");
		fp2 = fopen(argv[2],"r");
	}
	else{
		cerr << "default file" << endl;
		fp1 = fopen("./fighter/32_boundary.txt","r");
		fp2 = fopen("./fighter/32_corners.txt","r");
	}
	if(fp1 == NULL || fp2 == NULL){
		cerr << "fopen failed" << endl;
		return -1;
	}
	Points pt;
	pt.readPoints(fp1);
	pt.readCornerIdx(fp2);
	fclose(fp1);
	fclose(fp2);

	
	// 最適化
	cerr << "Points " << endl;
	Points pl = pt.rearrenge();
//	pl.print();
	cerr << "G" << endl;
	Graph G(pl);
	cerr << "C " << pt.corner.size() << endl;
	cerr << "Path" << endl;
	Points res = G.GetAllPoints();

	// スクリーンサイズに変形 (-1<x,y<1)
	RT r1 = pl.GetRect();
	RT r2 = pl.GetRect();
	RT rect;
	rect.top = max(r1.top,r2.top);
	rect.bottom = min(r1.bottom,r2.bottom);
	rect.left = min(r1.left,r2.left);
	rect.right = max(r1.right,r2.right);
	double range = max(rect.top-rect.bottom,rect.right-rect.left);
	Vec2 rectCenter((rect.right+rect.left)/2.0,(rect.top+rect.left)/2.0);
	pl.reshape(range,rectCenter);
	res.reshape(range,rectCenter);

	//出力
	double c1[] = { 1.0,0.0,0.0,1.0 };
	double c2[] = { 0.0,0.0,1.0,1.0 };
	cout << "fits = [" << endl;
	for(int i = 0; i < res.size(); i++){
		cout << res.at(i).x << ", " << res.at(i).y << ", " << "0.0," << endl;
	}
	cout << "];" << endl;
	cout << "color_f = [" << endl;
	for(int i = 0; i < res.size(); i++){
		cout << c1[0] << ", " << c1[1] << ", " << c1[2] << ", " << c1[3] << "," << endl;
	}
	cout << "];" << endl;
	cout << "pres = [" << endl;
	for(int i = 0; i < pl.size(); i++){
		cout << pl.at(i).x << ", " << pl.at(i).y << ", " << "0.0," << endl;
	}
	cout << "];" << endl;
	cout << "color_p = [" << endl;
	for(int i = 0; i < pl.size(); i++){
		cout << c2[0] << ", " << c2[1] << ", " << c2[2] << ", " << c2[3] << "," << endl;
	}
	cout << "];" << endl;
	/*
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
*/
	return 0;
}
