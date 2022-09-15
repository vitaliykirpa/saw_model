
CC=gcc
CXX=g++
RM=rm -f


CPPFLAGS=-std=c++14 -O3
LDFLAGS=-std=c++14 -O3
LDLIBS=

PROGRAM=mcmc

SRCS=Lattice.cpp \
     Model.cpp \
     Monte_Carlo.cpp


OBJS=$(subst .cpp,.o,$(SRCS))


all: $(PROGRAM)

mcmc.o: Lattice.cpp Model.cpp

main.o: Monte_Carlo.cpp


$(PROGRAM): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(LDLIBS)

clean:
	rm -rf $(OBJS) $(PROGRAM)
