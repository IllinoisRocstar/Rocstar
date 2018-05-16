# !/bin/csh

#
# Compiler to Use
#
FFF       = f90

#
# Any flags? ... comment out ones you don't want
#  -B101         # Generate code that will perform best on a Pentium Pro but still execute on a Pentium.
#  -P            # For profiling.  Run 'gprof fem_code > time'  to get a 'time' file
#  -g            # For Debugging.  Run 'gdb fem_code' and then type 'run' to debug
#  -O            # For Optimization
#
#TURING
FFLAGS     =  -lU77 -B101 -O

# Origin2000
# FFLAGS = -O3

#
# Object files for scale
#


OBJ     = data_declarations.o\
	    READ_FILES.o\
	    FIND_NODES.o\
	    TRANSLATE.o\
	    WRITE_OUTPUT.o\
	    WRITE_OUTPUT_2.o\
	    READ_IO_FILES.o\
	    READ_FRAC.o\
	    scale.o

#
#  Main executable files.
#
scale : $(OBJ)
	$(FFF) $(FFLAGS) $(OBJ) -o scale

#
# Dependancies
#
scale.o: scale.f90
	$(FFF) $(FFLAGS) scale.f90 -c -o scale.o

data_declarations.o: data_declarations.f90
	$(FFF) $(FFLAGS) data_declarations.f90 -c -o data_declarations.o	

READ_FILES.o: READ_FILES.f90
	$(FFF) $(FFLAGS) READ_FILES.f90 -c -o READ_FILES.o

FIND_NODES.o: FIND_NODES.f90
	$(FFF) $(FFLAGS) FIND_NODES.f90 -c -o FIND_NODES.o

TRANSLATE.o: TRANSLATE.f90
	$(FFF) $(FFLAGS) TRANSLATE.f90 -c -o TRANSLATE.o

WRITE_OUTPUT.o: WRITE_OUTPUT.f90
	$(FFF) $(FFLAGS) WRITE_OUTPUT.f90 -c -o WRITE_OUTPUT.o

READ_IO_FILES.o: READ_IO_FILES.f90
	$(FFF) $(FFLAGS) READ_IO_FILES.f90 -c -o READ_IO_FILES.o

WRITE_OUTPUT_2.o: WRITE_OUTPUT_2.f90
	$(FFF) $(FFLAGS) WRITE_OUTPUT_2.f90 -c -o WRITE_OUTPUT_2.o

READ_FRAC.o: READ_FRAC.f90
	$(FFF) $(FFLAGS) READ_FRAC.f90 -c -o READ_FRAC.o

#
#  Clean out directory by removing all object and executable files.
#
clean:
	/bin/rm *.o *~ *.mod scale

cleanobj:
	/bin/rm *.o

cleanmod:
	/bin/rm *.mod

cleanexe:
	/bin/rm scale

#cleantil:
#v	/bin/rm *~

