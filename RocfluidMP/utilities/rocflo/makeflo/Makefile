# Makefile for CSAR Rocflo partitioner
# Orion Sky Lawlor, olawlor@acm.org, 2003/7/17

# Set COMHOME only if not yet set, so that the makefile can be included 
# by the makefile in the Rocflo/util/makeflo directory
CFINTER = yes
COMHOME?=../../../../Roccom
LIBDIR?=$(COMHOME)/../lib
include $(COMHOME)/Makefile.common

##################################
# System-specific configuration

#OPTS            = -O

##################################
# Program objects, headers, and executables:

READOBJ=read.o readgrd.o readreg.o readmsh.o readhdf.o
WRITEOBJ=write.o writemsh.o writehdf.o writegrda.o writegrdb.o

UTIL_OBJ=hash_table.o read.a write.a util.o
OBJ= adj.o face.o buildface.o split.o volume.o \
	readinp.o writeflo.o writetop.o writemblock.o \
	$(UTIL_OBJ)
HEADERS=adj.h face.h util.h hash_table.h vector3d.h gridutil.h
DEST=$(BINPRE)makeflo    # Prepend $(BINPRE) to put into bin directory
RTST= $(UTIL_OBJ)
ATST=read.o $(UTIL_OBJ)

CPPFLAGS+=-DUSE_HDF=1 $(HDFINC)
LFLAGS=$(LDFLAGS) $(HDFLIB) $(SYSLIBS) $(DEPLIBS)

########################################

all:  $(DEST)

$(BINPRE)makeflo: makeflo.o  $(OBJ) $(HEADERS)
	$(LD) $(OPTS) $< $(OBJ) -o $@ $(LFLAGS)

$(BINPRE)makemblock: makemblock.o $(OBJ) $(HEADERS) 
	$(LD) $(OPTS) $< $(OBJ) -o $@ $(LFLAGS)

$(BINPRE)simpleblock: simpleblock.o $(OBJ)
	$(LD) $(OPTS) $< $(OBJ) -o $@ $(LFLAGS)

read.a: $(READOBJ)
	$(AR) $@ $(READOBJ)
	$(RANLIB) $@

write.a: $(WRITEOBJ)
	$(AR) $@ $(WRITEOBJ)
	$(RANLIB) $@

.SUFFIXES: .f90 .cpp

# Unit test programs
read.tst: $(RTST) read.cpp
	$(CXX)  $(OPTS) $(CCFLAGS) -DSTANDALONE -c read.cpp 
	$(LD) $(OPTS) read.o $(RTST) -o $@ $(LFLAGS)

adj.tst: $(ATST) adj.cpp
	$(CXX)  $(OPTS) $(CCFLAGS) -DSTANDALONE -c adj.cpp 
	$(LD) $(OPTS) adj.o $(ATST) -o $@ $(LFLAGS)

clean:
	-/bin/rm -fr *.o *.a *~ *.tst $(BINPRE)makeflo $(BINPRE)makemblock $(BINPRE)simpleblock ii_files/ SunWS_cache/
