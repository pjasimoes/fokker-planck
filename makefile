# makefile for fp in C++ program (PJSimoes, 2011)

CC=g++

CFLAGS=-Wall -c -O3
CFLAGS2=-Wall -c -pg
COPENMP=-fopenmp
CLIBS=-pg
EXECUTABLE=-o fp
SOURCES=fp.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all:	
	$(CC) $(CFLAGS) $(SOURCES)
	$(CC) $(CLIBS) $(OBJECTS) $(EXECUTABLE)

openmp:	
	$(CC) $(CFLAGS) $(COPENMP) $(SOURCES)
	$(CC) $(CLIBS) $(COPENMP) $(OBJECTS) $(EXECUTABLE)

single:	
	$(CC) $(CFLAGS) $(SOURCES)
	$(CC) $(CLIBS) $(OBJECTS) $(EXECUTABLE)

debug:	
	$(CC) $(CFLAGS2) $(SOURCES)
	$(CC) $(CLIBS) $(OBJECTS) $(EXECUTABLE)

clean:
		rm -rf *.o

#
#g++ -Wall -I/usr/local/include drive.cpp -c	
#g++ drive.o -lgsl -lgslcblas -lm -o drive
