CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3 -g -I. -DDASHMM_TESTING
INCLUDE = $(shell pkg-config --cflags hpx)
LIBS = $(shell pkg-config --libs hpx)

SRC = $(shell ls src/*.cc) testmain.cc
OBJ = $(SRC:.cc=.o)
EXEC = testmain

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) -o $(EXEC) $(OBJ) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDE)

clean:
	rm -rf $(EXEC) $(OBJ) *~
