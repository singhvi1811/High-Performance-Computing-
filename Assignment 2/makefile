SOURCES = precision.f90 lsq.f90 main.f90
OBJECTS = precision.o lsq.o main.o
MODULES = precision.mod lsq.mod
FC = ifort
OPT = 
LIBS = -mkl
PARALLEL = -qopenmp -heap-arrays 100  -traceback
DEBUG = -g
OTHER = -c
FFLAGS = $(OPT) $(OTHER) $(PARALLEL)

%.o: %.f90
	$(FC) $(FFLAGS) $<
%.o: %.f
	$(FC) $(FFLAGS) $<

all:  test_henon

test_henon: $(OBJECTS)
	$(FC) $(PARALLEL) $(OBJECTS) $(LIBS) -o test_henon 
	 
	


clean:
	rm -rf $(OBJECTS) $(MODULES) a.out core

