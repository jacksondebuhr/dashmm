CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3 -g
INCLUDE = -I../../. $(shell pkg-config --cflags hpx)
LIBS = ../../libdashmm.a $(shell pkg-config --libs hpx)

SRC = usermain.cc
OBJ = $(SRC:.cc=.o)

EXEC = userdemo

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) -o $(EXEC) $(OBJ) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDE)

clean:
	rm -rf $(EXEC) $(OBJ) *~
