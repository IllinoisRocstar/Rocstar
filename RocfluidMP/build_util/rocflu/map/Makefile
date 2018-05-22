################################################################################
#
# $Id: Makefile,v 1.6 2008/11/26 15:15:11 mtcampbe Exp $
#
# Purpose: Makefile for ROCFLU-MP's preprocessor.
#
# Description: none
#
# Notes: none
#
# Copyright: (c) 2001 by the University of Illinois
#
################################################################################

FLUIDHOME=../../..
SRCDIR = $(FLUIDHOME)/utilities/rocflu/map

EXEC=$(FLUIDHOME)/rflumap

LIBINT=$(FLUIDHOME)/$(BUILDLIB_DIR)/libflu.a

MODDIR=$(FLUIDHOME)/$(BUILDLIB_DIR)
SRCF77=
SRCF90=
SYSTEM=$(shell uname -s)

### Include common and dependent makefiles #####################################
include $(FLUIDHOME)/Makefile.common
include $(FLUIDHOME)/Makefile.$(SYSTEM)
LIBEXT= $(LIBSYS)

### Define SRCF90 and SRCF77 ###################################################
include $(SRCDIR)/$(FILELIST)

vpath %.F90 $(SRCDIR)
vpath Mod%.o .:$(FLUIDHOME)/$(BUILDLIB_DIR)
vpath RFLU_Mod%.o .:$(FLUIDHOME)/$(BUILDLIB_DIR)
vpath %.C $(SRCDIR)

include $(FLUIDHOME)/Makefile.dep

### Rules ######################################################################
OBJF77=$(SRCF77:.f=.o)
OBJF90=$(SRCF90:.F90=.o)
OBJCXX=$(SRCCXX:.C=.o)

### Targets ####################################################################
.PHONY:	all clean

all:    $(EXEC)

# Automatically determine the Fortran name mangling
main.o: rflumap.o
UNDERSCORE = $(if $(shell nm rflumap.o | \
             grep -i rflumap_),-DFORTRAN_UNDERSCORE,)
LOWERCASE =  $(if $(shell nm rflumap.o | \
             grep rflumap),-DFORTRAN_LOWERCASE,)
CXXFLAGS  += $(LOWERCASE) $(UNDERSCORE)

$(EXEC): $(OBJF77) $(OBJF90) $(LIBINT) $(OBJCXX)                    
	$(LINKER) $(OBJF77) $(OBJF90) $(OBJCXX) -o $(EXEC) $(LDFLAGS) $(LIBINT) $(LIBEXT) $(CXXLINKFLAGS)

clean:  
	$(RM) $(EXEC)
	$(RM) $(MODULES)
	$(RM) $(RIFFRAFF)

install:
	-$(CP) $(EXEC) $(BINDIR)/.

################################################################################
#
# RCS Revision history:
#
#   $Log: Makefile,v $
#   Revision 1.6  2008/11/26 15:15:11  mtcampbe
#   Port to NCSA Abe cluster for Rocflu.
#
#   Revision 1.5  2005/05/03 03:04:03  haselbac
#   Changed to linking with C++
#
#   Revision 1.4  2003/12/04 03:19:28  haselbac
#   Modified clean target
#
#   Revision 1.3  2003/03/31 19:51:44  haselbac
#   Removed Charm version of executable
#
#   Revision 1.2  2003/03/21 23:07:39  haselbac
#   Added LIBINT to deps, made CP robust, added _charm to EXEC
#
#   Revision 1.1.1.1  2003/03/20 22:48:04  haselbac
#   Initial revision
#
#   Revision 1.31  2003/02/02 21:17:11  haselbac
#   Inelegant bug fix for missing dependencies
#
#   Revision 1.30  2003/02/01 02:46:06  jiao
#   Fixed DEPPRE definition
#
#   Revision 1.29  2003/02/01 00:31:56  haselbac
#   Added RFLU_USER_EnforcePatchCoords
#
#   Revision 1.28  2003/01/31 22:12:07  jiao
#   Changed to compile modules in modfloflu to avoid duplication of object files.
#
#   Revision 1.27  2003/01/28 18:00:28  haselbac
#   Added and removed some files
#
#   Revision 1.26  2003/01/09 20:50:55  haselbac
#   Added RM variable
#
#   Revision 1.25  2002/12/27 22:07:14  jblazek
#   Splitted up RFLO_ModInterfaces and ModInterfaces.
#
#   Revision 1.24  2002/10/27 19:21:46  haselbac
#   Added RFLU_WriteGridInterface
#
#   Revision 1.23  2002/10/21 14:31:17  haselbac
#   Added rfluprep_charm binary if compiled with Charm
#
#   Revision 1.22  2002/10/16 21:21:54  haselbac
#   Added init for CHARM
#
#   Revision 1.21  2002/10/12 21:31:01  haselbac
#   Changed CHARM to RFLU_CHARM
#
#   Revision 1.20  2002/10/05 19:45:14  haselbac
#   Resolved conflict - CVS got confused
#
#   Revision 1.19  2002/09/27 00:57:11  jblazek
#   Changed makefiles - no makelinks needed.
#
#   Revision 1.18  2002/09/25 17:27:38  jferry
#   changed LIBS to SYSLIB
#
#   Revision 1.17  2002/09/17 22:51:24  jferry
#   Removed Fast Eulerian particle type
#
#   Revision 1.16  2002/09/09 16:37:37  haselbac
#   Added RFLU_ModConvertCv
#
#   Revision 1.15  2002/09/06 20:33:27  f-najjar
#   Move ModPartLag before ModBndPatch for Tile datastructure
#
#   Revision 1.14  2002/08/07 21:07:55  jblazek
#   Executables are left in the top-level directory instead of $(HOME)/bin.
#
#   Revision 1.13  2002/07/25 15:22:48  haselbac
#   Added RFLU_ModFEM
#
#   Revision 1.12  2002/06/27 15:57:00  haselbac
#   Added RFLU_ModCellFaceInfo and RFLU_ModCellMapping
#
#   Revision 1.11  2002/06/14 20:28:33  haselbac
#   Added ModSortSearch, RFLU_ModFaceList, RFLU_ModHashTable, deleted ModLocal
#
#   Revision 1.10  2002/06/10 21:37:43  haselbac
#   Added entries for BuildVersionString and RFLU_PrintHeader
#
#   Revision 1.9  2002/05/28 14:04:09  haselbac
#   For CHARM, add RFLU_CreatePartitions to SRC
#
#   Revision 1.8  2002/05/04 17:15:23  haselbac
#   Added ModLocal
#
#   Revision 1.7  2002/04/11 19:05:31  haselbac
#   Added new routines
#
#   Revision 1.6  2002/03/27 15:53:59  haselbac
#   Added entries for reading of input files
#
#   Revision 1.5  2002/03/01 17:10:00  haselbac
#   Deleted entry for LIBFLU_ModInterfaces
#
#   Revision 1.4  2002/02/08 14:24:07  haselbac
#   Added new routines for reading grids
#
#   Revision 1.3  2002/01/14 20:20:05  haselbac
#   Updates for new routines
#
#   Revision 1.2  2001/12/04 16:43:28  jblazek
#   Makefiles modified because the modules Global, BndPatch and Grid moved
#   to the modfloflu directory.
#
#   Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
#   Import of RocfluidMP
#
################################################################################
