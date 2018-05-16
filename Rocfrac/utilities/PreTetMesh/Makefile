SHELL	= /bin/sh

OSTYPE	= $(shell uname)

F90 = f90

OBJS = meshdata.o read_deck.o pat2tetmesh.o read_patran.o TetmeshOut.o

##if "$(OSTYPE)" "IRIX64"
#FFLAGS  = -64 -mips4 -extend_source -O3 -OPT:Olimit=0 
#BASICLIB = 
#LIBPATH = 
##endif

##ifeq "$(OSTYPE)" "Linux"
FFLAGS  = -W132 -O
BASICLIB = 
LIBPATH = 
##endif

.f.o:
	$(F90) $(FFLAGS) -c $*.f

pat2tetmesh: $(OBJS) 
	$(F90) $(FFLAGS) -o $@ $(OBJS) $(LIBPATH) $(BASICLIB)

clean:
	rm -f *.o pat2tetmesh *.mod 