all : test 

test : test.cpp
	g++ $^ -std=c++11 -O2 -o test;

