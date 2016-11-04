CXX = mpicxx
CXXFLAGS = -std=c++11 -Wall -O3 -g
INCLUDE = $(shell pkg-config --cflags hpx)
LIBS = $(shell pkg-config --libs hpx)

SRC = $(shell ls *.cc)
OBJ = $(SRC:.cc=.o)
EXEC = treehpxv2

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) -o $(EXEC) $(OBJ) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDE)

clean:
	rm -rf $(EXEC) $(OBJ) *~
