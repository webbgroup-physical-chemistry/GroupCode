CC=g++
CFLAGS=-I. -O
DEPS = MyOptionParser.h prob.h
OBJ = main.o MyOptionParser.o prob.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(OBJ)
	$(CC) -o Boltzmann_Weight $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(OBJ) Boltzmann_Weight
