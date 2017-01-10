CXX = g++
CXXFLAGS = -std=c++11 -Wall -g -O0

INCLUDE = 
LIBS = 

SRC = $(shell ls *.cc)
OBJ = $(SRC:.cc=.o)
EXEC = convtoadjmat

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) -o $(EXEC) $(OBJ) $(LIBS)

%.o:%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDE)

clean:
	rm -rf $(EXEC) $(OBJ) *~
