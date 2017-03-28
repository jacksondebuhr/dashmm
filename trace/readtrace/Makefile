CXX = g++
CXXFLAGS = -std=c++11 -Wall -O0 -g
INCLUDE = -I/usr/local/include `freetype-config --cflags`
LIBS = -L/usr/local/lib -lpng -lpngwriter -lz -lfreetype -lsqlite3

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
