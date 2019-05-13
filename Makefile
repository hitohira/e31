all : test 

test : test.cpp primitives.h graph.h
	g++ $< -std=c++11 -O2 -o test;

