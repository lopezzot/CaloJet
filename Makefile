ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) 
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      =-I$(ROOTSYS)/include -O -Wall -fPIC
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

FASTJET_CONFIG=/home/software/fastjet/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
CXXFLAGS    += $(ROOTCFLAGS)
CXXFLAGS    += $(shell $(FASTJET_CONFIG) --cxxflags)
LIBS        = $(ROOTLIBS) $(LIBSFASTJET)
GLIBS       = $(ROOTGLIBS)

all: dohist

dohist: dohist.o MyTree.o
	$(CXX) -o $@ dohist.o MyTree.o $(CXXFLAGS) $(LIBS)

dohist.o:  dohist.C 
	$(CXX) -c $(CXXFLAGS) -I.  -o dohist.o dohist.C

MyTree.o:  MyTree.C MyTree.h
	$(CXX) -c $(CXXFLAGS) -I.   -o MyTree.o MyTree.C

# clean
clean:
	rm -f *~ *.o *.o~ core

