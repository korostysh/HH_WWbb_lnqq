ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

#ROOTCFLAGS    = $(shell /usr/bin/root-config --cflags)
#ROOTLIBS      = $(shell /usr/bin/root-config --libs)
#ROOTGLIBS     = $(shell /usr/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -Wno-deprecated

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit

CXXFLAGS      += $(ROOTCFLAGS)
CXX           += -I./	
LIBS           = $(ROOTLIBS) 

GLIBS          = $(filter-out -lNew, $(NGLIBS))

CXX	      += -I./obj/
OUTLIB	      = ./obj/
.SUFFIXES: .C
.PREFIXES: ./obj/

#----------------------------------------------------#

lab1: lab1.cpp
	g++ -o lab1 lab1.cpp `root-config --cflags --glibs ` -lSpectrum

Plot: widthplot.cpp
	g++ -o widthplot widthplot.cpp `root-config --cflags --glibs ` -lSpectrum

convertL: convertL.cc
	g++ -o convertL convertL.cc `root-config --cflags --glibs`

plotsFast: plotsFast.cc analysis_corto.o
	g++ -o plotsFast plotsFast.cc *.o `root-config --cflags --glibs`

analysis_corto.o: analysis_corto.cc analysis_corto.h
	$(CXX) $(CXXFLAGS) -c -I. -o analysis_corto.o $<

clean:
	rm -f lab1
	rm -f Data
	rm -f widthplot
	rm -f convertL	rm -f 
	rm -f *~
	rm -f plots/*.png
	rm -f HME/Spread/*.png
	rm -f HME/No_spread/*.png
	rm -f HME/*.png
	rm -f HME/HME_both/*.png
	rm -f HME/HME2D/*.png
	rm -f HME/Det/*.png
	rm -f HME/CDF/*.png
	

cleananalysis:
	rm -f analysis_corto.o
	rm -f plotsFast
	rm -f *~
