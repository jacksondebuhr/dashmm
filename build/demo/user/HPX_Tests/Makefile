CXX = g++
CXXFLAGS = -std=gnu++11 -Wall -O3 -g
INCLUDE = $(shell pkg-config --cflags hpx)
LIBS = $(shell pkg-config --libs hpx)

SRC = HPX_Test.cc
OBJ = $(SRC:.cc=.o)

EXEC = HPX_Test

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) -o $(EXEC) $(OBJ) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDE)

clean:
	rm -rf $(EXEC) $(OBJ) *~
