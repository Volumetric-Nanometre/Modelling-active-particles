IDIR =../include
CC=gcc
CFLAGS=-Wall -g -I$(IDIR)

ODIR=obj
LDIR =../lib

LIBS=-lm

_DEPS = diffusionmatrix.h maths_functions.h particles.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = diffusionmatrix.o main.o maths_functions.o particles.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

../bin/executable: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~