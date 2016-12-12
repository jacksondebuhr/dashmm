CXX = mpicxx
CXXFLAGS = -std=c++11 -Wall -O0 -g
INCLUDE =
LIBS =

SRC = $(shell ls *.cc)
OBJ = $(SRC:.cc=.o)

EXEC = test

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) -o $(EXEC) $(OBJ) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDE)

clean:
	rm -rf $(EXEC) $(OBJ) *~
