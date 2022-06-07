IDIR =include
CC=g++
CFLAGS=-I$(IDIR) -Iinclude* -Wall -std=c++17 # -Wextra 

ODIR=obj

LIBS=-lm

_DEPS = SSG.h Strategy.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o SSG.o Strategy.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -g -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) -g -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
