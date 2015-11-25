mark:	lepto.cxx epsilon.h epsilon.c washout.c washout.h dispersion.c dispersion.h
	g++ -g -std=c++0x -c epsilon.c -o epsilon.o -I. -I/mt/home/mark/etc/include
	g++ -g -std=c++0x -c washout.c -o washout.o -I. -I/mt/home/mark/etc/include
	g++ -g -std=c++0x -c dispersion.c -o dispersion.o -I. -I/mt/home/mark/etc/include
	g++ -g -std=c++0x -c -o lepto.o  lepto.cxx -I. -I/mt/home/mark/etc/include
	g++ -g -std=c++0x -o lepto dispersion.o washout.o epsilon.o lepto.o -lgomp -lnlopt -lgsl -lgslcblas -L/usr/lib -L/mt/home/mark/etc/lib
jess:	lepto.cxx epsilon.h epsilon.c washout.c washout.h dispersion.c dispersion.h
	g++ -g -std=c++0x -c epsilon.c -o epsilon.o -I. -I/opt/local/include
	g++ -g -std=c++0x -c washout.c -o washout.o -I. -I/opt/local/include
	g++ -g -std=c++0x -c dispersion.c -o dispersion.o -I. -I/opt/local/include
	g++ -g -std=c++0x -c -o lepto.o  lepto.cxx -I. -I/opt/local/include
	g++ -g -std=c++0x -o lepto dispersion.o washout.o epsilon.o lepto.o -lgsl -lgslcblas -L/usr/lib -L/opt/local/lib

