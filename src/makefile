IDIR =../include
CC=gcc
CFLAGS=-Wall -std=c11 -fopenmp -O3 -I$(IDIR) -I"C:\Program Files (x86)\Microsoft SDKs\MPI\Include" -I"C:\Users\mieha\Documents\GitHub\Discrete-dipole-approximation\include"


ODIR=obj
LDIR =../lib

LIBS=-lm -lmosh-dda -lgsl -lgslcblas -lmsmpi  -lfftw3 -L$(LDIR) -L"C:\Program Files (x86)\Microsoft SDKs\MPI\Lib" -L"C:\Users\mieha\Documents\GitHub\Discrete-dipole-approximation\lib"

_DEPS = diffusionmatrix.h maths_functions.h particles.h stochastic_force.h moving_on.h forces.h initial_finalisation.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o diffusionmatrix.o maths_functions.o particles.o stochastic_force.o moving_on.o forces.o initial_finalisation.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

args =

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) ${args}

../bin/executable: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) ${args}

.PHONY: clean run

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

run:
	mpiexec -n $(n) ../bin/executable.exe ${args}
