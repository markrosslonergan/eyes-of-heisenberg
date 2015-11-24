d1:	lepto.cxx epsilon.h epsilon.c
	g++ -g -std=c++0x -c epsilon.c -o epsilon.o -I. -I/mt/home/mark/etc/include
	g++ -g -std=c++0x -c washout.c -o washout.o -I. -I/mt/home/mark/etc/include
	g++ -g -std=c++0x -c dispersion.c -o dispersion.o -I. -I/mt/home/mark/etc/include
	g++ -g -std=c++0x -c -o lepto.o  lepto.cxx -I. -I/mt/home/mark/etc/include
	g++ -g -std=c++0x -o lepto dispersion.o washout.o epsilon.o lepto.o -lgomp -lnlopt -lgsl -lgslcblas -L/usr/lib -L/mt/home/mark/etc/lib

