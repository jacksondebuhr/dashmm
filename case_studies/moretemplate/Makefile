CXX = g++
CXXFLAGS = -std=c++11 -Wall -O0 -g
INCLUDE = #$(shell pkg-config --cflags hpx)
LIBS = #$(shell pkg-config --libs hpx)

SRC = main.cc
OBJ = $(SRC:.cc=.o)

EXEC = template

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) -o $(EXEC) $(OBJ) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDE)

clean:
	rm -rf $(EXEC) $(OBJ) *~
