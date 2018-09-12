IDIR =../include
CC=gcc
CFLAGS=-Wall -fopenmp -O3 -I$(IDIR) -IC:/GSL/include


ODIR=../src/obj
LDIR =../lib

LIBS=-lm -lgsl -lgslcblas -lmsmpi -L$(LDIR) -LC:/GSL/lib

_DEPS = diffusionmatrix.h maths_functions.h particles.h stochastic_force.h moving_on.h forces.h initial_finalisation.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ =  test.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

args =

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) ${args}

../bin/unit-test: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) ${args}

.PHONY: clean run

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

run:
	mpiexec -n $(n) ../bin/unit-test.exe ${args}
