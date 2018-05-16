# Origin 2000

#F90 = f90
#FFLAGS     =  -extend_source  -O3 -OPT:Olimit=0

# Turing

F90 = f90
FFLAGS = -W132 -O -lU77

# Sun

#F90 = f90
#FFLAGS = -O -e


#
# Object files for cohesiv2
#
OBJF90 = CharKeyWrd.o


OBJ     = meshdata.o \
	linked_list.o \
	linked_list_2d.o \
	meshtran3d.v1.o \
	readinp.o \
	read_ansys.o \
	read_patran.o \
	read_patran_cohin.o \
	read_tetmesh_cohin.o \
	read_tetmesh.o \
	roc_intrfc_data.comm.o \
	STRUCT_tri.o\
	coh_opposite_side.o

all: ./MeshTran

#../meshtran3d.v1: $(OBJ)
#	$(F90) $(FFLAGS) $(OBJ) -o ../meshtran3d.v1 -lmetis -lU77 -L .


./MeshTran: $(OBJ) $(OBJF90)
	$(F90) $(FFLAGS) $(OBJ) $(OBJF90) -o ./MeshTran -L. -lmetis -lm

linked_list.o: linked_list.f
	$(F90) $(FFLAGS) linked_list.f -c -o linked_list.o

linked_list_2d.o: linked_list_2d.f
	$(F90) $(FFLAGS) linked_list_2d.f -c -o linked_list_2d.o

meshdata.o: meshdata.f
	$(F90) $(FFLAGS) meshdata.f -c -o meshdata.o

meshtran3d.v1.o: meshtran3d.v1.f meshdata.f
	$(F90) $(FFLAGS) meshtran3d.v1.f -c -o meshtran3d.v1.o

roc_intrfc_data.comm.o: roc_intrfc_data.comm.f meshdata.f
	$(F90) $(FFLAGS) roc_intrfc_data.comm.f -c -o roc_intrfc_data.comm.o

readinp.o: readinp.f
	$(F90) $(FFLAGS) readinp.f -c -o readinp.o

read_ansys.o: read_ansys.f meshdata.f
	$(F90) $(FFLAGS) read_ansys.f -c -o read_ansys.o

read_patran.o: read_patran.f meshdata.f
	$(F90) $(FFLAGS) read_patran.f -c -o read_patran.o

read_patran_cohin.o: read_patran_cohin.f meshdata.f
	$(F90) $(FFLAGS) read_patran_cohin.f -c -o read_patran_cohin.o

read_tetmesh_cohin.o: read_tetmesh_cohin.f meshdata.f
	$(F90) $(FFLAGS) read_tetmesh_cohin.f -c -o read_tetmesh_cohin.o

read_tetmesh.o: read_tetmesh.f meshdata.f
	$(F90) $(FFLAGS) read_tetmesh.f -c -o read_tetmesh.o

STRUCT_tri.o: STRUCT_tri.f
	$(F90) $(FFLAGS) STRUCT_tri.f -c -o STRUCT_tri.o

coh_opposite_side.o: coh_opposite_side.f
	$(F90) $(FFLAGS) coh_opposite_side.f -c -o coh_opposite_side.o

CharKeyWrd.o: ../../Source/CharKeyWrd.f90
	$(F90) -O -lU77  ../../Source/CharKeyWrd.f90 -c -o CharKeyWrd.o

clean:
	rm -f *.o *mod ./meshtran3d.v1



