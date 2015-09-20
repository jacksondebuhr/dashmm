CXX = icpc #g++

CXXFLAGS = -std=c++11 -Wall -pedantic -g -O3

INCLUDE = #$(shell pkg-config --cflags gsl)
LIBS = #$(shell pkg-config --libs gsl)

OPTIONS = #-pg

SRC = $(shell ls *.cc)
OBJ = $(SRC:.cc=.o)
EXEC = mesodashmm

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(OPTIONS) -o $(EXEC) $(OBJ) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(OPTIONS) -c $< -o $@ $(INCLUDE)

clean:
	rm -rf $(EXEC) $(OBJ) *~
