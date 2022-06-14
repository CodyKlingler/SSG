IDIR =include
CC=g++
CFLAGS=-I$(IDIR) -Iinclude* -Wall -std=c++17 -O0# -Wextra 

ODIR=obj

LIBS=-lm

_DEPS = SSG.h Strategy.h SSG_tests.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o SSG.o Strategy.o SSG_tests.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -g -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) -g -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
