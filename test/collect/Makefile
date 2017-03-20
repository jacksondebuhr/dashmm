CXX = mpicxx
CXXFLAGS = -std=gnu++11 -Wall -O3 -g #-DDASHMM_INSTRUMENTATION
INCLUDE = -I../../include $(shell pkg-config --cflags hpx)
LIBS = -L../../lib -ldashmm $(shell pkg-config --libs hpx)

SRC = collect.cc
OBJ = $(SRC:.cc=.o)

EXEC = collect

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) -o $(EXEC) $(OBJ) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDE)

clean:
	rm -rf $(EXEC) $(OBJ) *~
