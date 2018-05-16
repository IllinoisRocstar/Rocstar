PROGRAM = CohIn
FPROGRAM = fort_prg
LIBRARY = libcomesh.a	

TARGETS =  $(PROGRAM) $(FPROGRAM) $(LIBRARY) 

OPTS = -O2  #-g

LIBDIR = -L.
INCDIR = 
LIB= -lcomesh
LIB_SRC = MVec.cpp Face.cpp QuadFace.cpp TriFace.cpp Element.cpp  \
	HexElement.cpp TetElement.cpp TCoElement.cpp QCoElement.cpp \
	Node.cpp Mesh.cpp IntList.cpp NodeList.cpp FaceList.cpp \
	FaceListList.cpp ElementList.cpp 
PROG_SRC = mainc.c
LIB_OBJ = MVec.o Face.o QuadFace.o TriFace.o Element.o  \
	HexElement.o TetElement.o TCoElement.o QCoElement.o Node.o Mesh.o \
	IntList.o NodeList.o FaceList.o FaceListList.o ElementList.o
PROG_OBJ = mainc.o
FORT_OBJ = mainf.o

default: $(LIBRARY) $(PROGRAM) 

library: $(LIBRARY)

.cpp.o  :
	g++ $(OPTS) -c $(INCDIR) $<

.cc.o  :
	g++ $(OPTS) -c $(INCDIR) $<

.c.o  :
	g++ $(OPTS) -c $(INCDIR) $<

.f90.o  :
	f90 $(OPTS) -c $(INCDIR) $<

$(PROGRAM) : $(PROG_OBJ) $(LIBRARY)
	g++  $(OPTS) -o $(PROGRAM)  $(PROG_OBJ)  $(LIBDIR) $(LDFLAGS) $(LIB)

$(LIBRARY) : $(LIB_OBJ) $(LIB_SRC)
	ar rv $(LIBRARY) $(LIB_OBJ)
	ranlib $(LIBRARY)

$(FPROGRAM) : $(FORT_OBJ) $(LIBRARY)
	f90 $(OPTS) -o $(FPROGRAM) $(FORT_OBJ) $(LIBRARY) -L/usr/local/lib/gcc-lib/sparc-sun-solaris2.6/2.8.1/ -L/usr/local/lib -lgcc -lobjc -lstdc++

depend: $(LIB_SRC) $(PROJ_SRC) 
	makedepend $(LIB_SRC) $(PROJ_SRC) $(INCDIR)
# DO NOT DELETE THIS LINE -- make depend depends on it.

MVec.o: MVec.hpp /usr/include/math.h general.hpp
Face.o: Face.hpp Node.hpp MVec.hpp /usr/include/math.h general.hpp
Face.o: Element.hpp
QuadFace.o: QuadFace.hpp Face.hpp Node.hpp MVec.hpp /usr/include/math.h
QuadFace.o: general.hpp Element.hpp QCoElement.hpp
TriFace.o: TriFace.hpp Face.hpp Node.hpp MVec.hpp /usr/include/math.h
TriFace.o: general.hpp Element.hpp TCoElement.hpp
Element.o: Element.hpp Face.hpp Node.hpp MVec.hpp /usr/include/math.h
Element.o: general.hpp HexElement.hpp TetElement.hpp TCoElement.hpp
Element.o: QCoElement.hpp Mesh.hpp /usr/include/assert.h
HexElement.o: HexElement.hpp Element.hpp Face.hpp Node.hpp MVec.hpp
HexElement.o: /usr/include/math.h general.hpp QuadFace.hpp Mesh.hpp
TetElement.o: TetElement.hpp Element.hpp Face.hpp Node.hpp MVec.hpp
TetElement.o: /usr/include/math.h general.hpp TriFace.hpp Mesh.hpp
TCoElement.o: TCoElement.hpp Element.hpp Face.hpp Node.hpp MVec.hpp
TCoElement.o: /usr/include/math.h general.hpp TriFace.hpp Mesh.hpp
QCoElement.o: QCoElement.hpp Element.hpp Face.hpp Node.hpp MVec.hpp
QCoElement.o: /usr/include/math.h general.hpp QuadFace.hpp Mesh.hpp
Node.o: Node.hpp MVec.hpp /usr/include/math.h general.hpp Face.hpp
Node.o: Element.hpp FaceList.hpp ElementList.hpp Mesh.hpp
Node.o: /usr/include/assert.h
Mesh.o: Mesh.hpp general.hpp Element.hpp Face.hpp Node.hpp MVec.hpp
Mesh.o: /usr/include/math.h IntList.hpp /usr/include/assert.h NodeList.hpp
Mesh.o: ElementList.hpp FaceList.hpp FaceListList.hpp
IntList.o: IntList.hpp general.hpp /usr/include/assert.h
NodeList.o: NodeList.hpp general.hpp Node.hpp MVec.hpp /usr/include/math.h
FaceList.o: FaceList.hpp general.hpp Face.hpp Node.hpp MVec.hpp
FaceList.o: /usr/include/math.h Element.hpp
FaceListList.o: FaceListList.hpp general.hpp FaceList.hpp Face.hpp Node.hpp
FaceListList.o: MVec.hpp /usr/include/math.h Element.hpp
ElementList.o: ElementList.hpp general.hpp Element.hpp Face.hpp Node.hpp
ElementList.o: MVec.hpp /usr/include/math.h
