CXX = mpicxx
CXXFLAGS = -std=c++11 -Wall -g -O0

INCLUDE = -I../../include $(shell pkg-config --cflags hpx)
LIBS = -I../../lib -ldashmm $(shell pkg-config --libs hpx)

SRC = $(shell ls *.cc)
OBJ = $(SRC:.cc=.o)
EXEC = test

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) -o $(EXEC) $(OBJ) $(LIBS)

%.o:%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDE)

clean:
	rm -rf $(EXEC) $(OBJ) *~
