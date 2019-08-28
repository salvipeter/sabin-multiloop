all: sabin

TRANSFINITE=/home/salvi/project/transfinite
INCLUDES=-I$(TRANSFINITE)/src/geom -I/usr/include/eigen3
CXXFLAGS=-g -Wall -std=c++17 $(INCLUDES)
LIBS=-L$(TRANSFINITE)/release/geom -lgeom

sabin: sabin.o solver.o
	g++ -o $@ $^ $(LIBS)
