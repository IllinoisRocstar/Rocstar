# Origin 2000

F90 = f90
FFLAGS     =  -extend_source  -O3 -OPT:Olimit=0

# Turing

#F90 = f90
#FFLAGS = -W132 -O 


#
# Object files 
#

OBJ	= meshdata.o \
	CohIn2Fract.o \
	read_deck.o

all: ./CohIn2Fract

CohIn2Fract:  $(OBJ)
	$(F90) $(FFLAGS) $(OBJ) -o CohIn2Fract

# Dependancies

meshdata.o: ../PreTetMesh/meshdata.f
	$(F90) $(FFLAGS) ../PreTetMesh/meshdata.f -c -o meshdata.o

CohIn2Fract.o: CohIn2Fract.f ../PreTetMesh/meshdata.f
	$(F90) $(FFLAGS) CohIn2Fract.f -c -o CohIn2Fract.o

read_deck.o:../PreTetMesh/read_deck.f ../PreTetMesh/meshdata.f
	$(F90) $(FFLAGS) ../PreTetMesh/read_deck.f -c -o read_deck.o