all : test 

test : test.cpp
	g++ $^ -std=c++11 -isystem ./eigen
