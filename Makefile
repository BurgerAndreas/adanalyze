#DEBUG ?= 0
#ifeq (DEBUG, 1)
CXXFLAGS= -pedantic -fopenmp -O0 -g  -Wall -std=c++11 
LDFLAGS= -g -L/usr/lib -lgsl -lgslcblas -fopenmp
#else
#CXXFLAGS= -std=c++11 -pedantic  -fopenmp -O3 -g -Wall 
#LDFLAGS= -fopenmp -g -lgsl -lgslcblas 
#endif

CC = g++

fsim.o: fsim.cpp vfield.hpp Makefile

vfield.o: vfield.cpp vfield.hpp Makefile

vfield_ww.o: vfield_ww.cpp vfield.hpp Makefile

vfield_dyn.o: vfield_dyn.cpp vfield.hpp vfield_ww.cpp Makefile

delabella.o: delabella.cpp delabella.h Makefile

ad_analyze.o: ad_analyze.cpp delabella.cpp delabella.h vfield.hpp vfield_ww.cpp Makefile

fsim: fsim.o vfield.o vfield_ww.o vfield_dyn.o ad_analyze.o delabella.o Makefile
	g++ -fopenmp -Wall -o fsim fsim.o vfield.o vfield_ww.o vfield_dyn.o ad_analyze.o delabella.o -lgsl -lgslcblas 

clean : *.o
	rm *.o fsim tfsim

fsim_act.tgz : fsim.cpp vfield.cpp vfield_ww.cpp vfield_dyn.cpp vfield.hpp Makefile
	tar -czf fsim_act.tgz fsim.cpp vfield.cpp vfield_ww.cpp vfield_dyn.cpp vfield.hpp Makefile

asim_act.tgz : fsim.cpp vfield.cpp vfield_ww.cpp vfield_dyn.cpp vfield.hpp Makefile analyze/Makefile analyze/avfield.hpp analyze/avfield.cpp analyze/asim.cpp 
	tar -czf asim_act.tgz fsim.cpp vfield.cpp vfield_ww.cpp vfield_dyn.cpp vfield.hpp Makefile analyze/Makefile analyze/avfield.hpp analyze/avfield.cpp analyze/asim.cpp 
