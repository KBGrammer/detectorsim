CCC         = mpicxx
OPT         = -O3

all: detectorsim

vector3.o: vector3.cpp vector3.h
	$(CCC) $(OPT) -c vector3.cpp

rng.o: rng.cpp rng.h
	$(CCC) $(OPT) -c rng.cpp

detectorsim: main.cpp vector3.o rng.o
	$(CCC) $(OPT) -o detectorsim main.cpp vector3.o rng.o

profile: detectorprof

vector3p.o: vector3.cpp vector3.h
	$(CCC) -pg -c vector3.cpp

rngp.o: rng.cpp rng.h
	$(CCC) -pg -c rng.cpp

detectorprof: main.cpp vector3p.o rngp.o
	$(CCC) -pg -o detectorsim main.cpp vector3.o rng.o

clean:
	rm -f *.o *~ */*.e4* */*.pe4* */*.po4*
