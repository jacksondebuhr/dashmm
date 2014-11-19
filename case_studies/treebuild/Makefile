#for frigg
CC       =  mpicc -std=gnu99

#for cutter
#CC       =  mpicc -std=gnu99

EXEC   = test-driver 

OBJS   = dashmm-tree.o \
	test-driver.o

INCL   = dashmm-tree.h \
	Makefile

#Some other compilation options
PROFILE = #-pg
DEBUGGING = -g


#Some HPX settings

#for frigg
HPXINCLUDE = -I/u/jdebuhr/Code/hpx/build/include
HPXLIBS = -L/u/jdebuhr/Code/hpx/build/lib -lhpx -lpthread

#for cutter
#HPXINCLUDE = -I/u/jdebuhr/hpx/build/include
#HPXLIBS = -L/u/jdebuhr/hpx/build/lib -lhpx -lrdmacm -libverbs

CFLAGS = -O3 -Wall -pedantic $(DEBUGGING) $(PROFILE) $(HPXINCLUDE)


LIBS   = -lm $(HPXLIBS)


$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS) -o  $(EXEC) $(PROFILE)

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC)


