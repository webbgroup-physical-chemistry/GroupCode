CC=g++
CFLAGS=-I. -O
DEPS = grid.h read_input.h read_pqr.h MyOptionParser.h
OBJ = main.o grid.o read_input.o read_pqr.o MyOptionParser.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(OBJ)
	$(CC) -o Dielectric_Map $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(OBJ) Dielectric_Map
