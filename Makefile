directory:
	mkdir -p bin
	mkdir -p obj

all:directory third
	g++ src/Dense.cpp -o obj/Dense.o -c -O3 -std=c++11 -flto -I ./third/ -I/usr/include/eigen3
	g++ obj/*.o -o bin/Dense `pkg-config --cflags --libs opencv` -flto

third:directory
	g++ third/GaussianPyramid.cpp -o obj/GaussianPyramid.o -c -O3 -std=c++11 -flto
	g++ third/Stochastic.cpp -o obj/Stochastic.o -c -O3 -std=c++11 -flto
	g++ third/OpticalFlow.cpp -o obj/OpticalFlow.o -c -O3 -std=c++11 -flto

clean:
	rm -f obj/*
	rm -f bin/*

