################################################################################
#
# $Id: Makefile,v 1.13 2008/11/26 15:15:09 mtcampbe Exp $
#
# Purpose: Makefile for RocfluidMP standalone codes
#
# Copyright: (c) 2003 by the University of Illinois
#
################################################################################

ifdef RFLO
  EXEC=../rflomp
  LIBINT= ../$(BUILDLIB_DIR)/libflo.a
endif
ifdef RFLU
  LIBINT=../$(BUILDLIB_DIR)/libflu.a
  EXEC=../rflump
endif
FLUIDHOME=../
MODDIR=../$(BUILDLIB_DIR)
SRCF77=
SRCF90=
SYSTEM=$(shell uname -s)

FPPFLAGS+= -I../modflo

### Include common and dependent makefiles #####################################
include ../Makefile.common
include ../Makefile.$(SYSTEM)
LIBEXT= $(LIBSYS) $(LIBMPI) $(LIBLAP) $(LIBPET) $(LIBX11)

### vpath definitions ##########################################################
vpath Mod%.F90 		../modfloflu

ifdef RFLO
  vpath RFLO_Mod%.F90	../modflo
  vpath RFLO_%.F90	../libflo
  vpath %.F90		../libfloflu:../rocflo:../standalone
endif
ifdef RFLU
  vpath RFLU_Mod%.F90	../modflu
  vpath RFLU_%.F90	../libflu
  vpath %.F90           ../libfloflu:../rocflu:../standalone
  vpath %.f		../libflu
  vpath %.C		../standalone
endif

vpath %.o .:../$(BUILDLIB_DIR)

### Define SRCF90 and SRCF77 ###################################################
include ../standalone/$(FILELIST)

include ../Makefile.dep

### Rules ######################################################################
OBJF77=$(SRCF77:.f=.o)
OBJF90=$(SRCF90:.F90=.o)
OBJCXX=$(SRCCXX:.C=.o)

### Targets ####################################################################
.PHONY:	all clean

all:    $(EXEC)

  ifdef RFLU
main.o: rflump.o
  # Automatically determine the Fortran name mangling
    UNDERSCORE = $(if $(shell nm rflump.o | \
                 grep -i rflump_),-DFORTRAN_UNDERSCORE,)
    LOWERCASE = $(if $(shell nm rflump.o | \
                grep rflump),-DFORTRAN_LOWERCASE,)
    CXXFLAGS  += $(LOWERCASE) $(UNDERSCORE)
  endif

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
# $Log: Makefile,v $
# Revision 1.13  2008/11/26 15:15:09  mtcampbe
# Port to NCSA Abe cluster for Rocflu.
#
# Revision 1.12  2005/07/05 16:15:18  haselbac
# Removed CHARM ifdef, cosmetics
#			        
# Revision 1.11  2005/05/03 03:03:21  haselbac  
# Changed to linking with C++		        
#
# Revision 1.10  2004/12/19 15:42:38  haselbac  
# Expanded LIBEXT list  		        
#
# Revision 1.9  2003/12/04 03:19:18  haselbac   
# Modified clean target 		        
#
# Revision 1.8  2003/11/20 16:40:32  mdbrandy   
# Backing out RocfluidMP changes from 11-17-03  
#
# Revision 1.5  2003/05/15 02:56:59  jblazek    
# Inlined index function.		        
#
# Revision 1.4  2003/04/11 15:15:18  jiao       
# Changed RFLU_CHARM back to CHARM.	        
#
# Revision 1.3  2003/03/21 23:05:34  haselbac   
# Added LIBINT to deps, made CP robust          
#
# Revision 1.2  2003/03/21 17:09:09  haselbac   
# Added missing header and footer	        
#
################################################################################
