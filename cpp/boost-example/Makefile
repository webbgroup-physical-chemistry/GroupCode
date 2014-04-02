NAME := integrate
EXTRA_SRC := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(EXTRA_SRC))

CC := /opt/apps/intel/13/composer_xe_2013.2.146/bin/intel64/icc 
INCLUDE := -I. -I$(TACC_GSL_INC) -I$(TACC_GSL_INC)/gsl -I$(TACC_BOOST_INC)
LDFLAGS := -L$(TACC_GSL_LIB) -L$(TACC_BOOST_LIB)
LDLIBS := -lgsl -lgslcblas -lboost_program_options 
CPPFLAGS := $(INCLUDE)

$(NAME): $(OBJS)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^
%.o: %.cpp
	$(CC) $(CPPFLAGS) $(INCLUDE) -c $<
clean:
	rm -f $(NAME) $(OBJS)

.PHONY: clean
