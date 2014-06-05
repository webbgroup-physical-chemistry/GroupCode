NAME := WHAM
EXTRA_SRC := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(EXTRA_SRC))

BOOST_INC := /usr/local/Cellar/boost/1.55.0_1/include
BOOST_LIB :=  /usr/local/Cellar/boost/1.55.0_1/lib

HDF5_INC := /usr/local/Cellar/hdf5/1.8.8/include
HDF5_LIB := /usr/local/Cellar/hdf5/1.8.8/lib
    
CC := g++
INCLUDE := -I. -I$(BOOST_INC) -I$(HDF5_INC)
LDFLAGS := -L$(BOOST_LIB) -L$(HDF5_LIB)
LDLIBS := -lboost_program_options -lhdf5_hl_cpp -lhdf5_cpp -lhdf5 -lz -lm
CPPFLAGS := $(INCLUDE) -O3 -Wall

$(NAME): $(OBJS)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^

%.o: %.cpp 
	$(CC) $(CPPFLAGS) -c $<

clean:
	rm -f $(NAME) $(OBJS)

.PHONY: clean
