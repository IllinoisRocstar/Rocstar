OBJ =	neighbor.o struct.o
#Origin 2k
#fc = f90 -extend_source
#opt = -O3

#Turing
#fc = f90 -W132
#opt = -O

#Sun
fc = f90 -e
opt = -O

tetnhbor:	$(OBJ)
	$(fc) $(opt) -o tetnhbor $(OBJ)

.f.o:	
	$(fc) $(opt) -c $*.f
clean:
	rm -f *.o tetnhbor
