IDIR =../include
CC=gcc
CFLAGS=-Wall -g -fopenmp -I$(IDIR)

ODIR=obj
LDIR =../lib

LIBS=-lm

_DEPS = diffusionmatrix.h maths_functions.h particles.h stochastic_force.h moving_on.h initial_velocities.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = diffusionmatrix.o main.o maths_functions.o particles.o stochastic_force.o moving_on.o initial_velocities.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

../bin/executable: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
