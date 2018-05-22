################################################################################
#
# $Id: Makefile,v 1.4 2008/11/26 15:15:10 mtcampbe Exp $
#
# Purpose: Makefile for ROCFLU-MP's initializor.
#
# Description: None
#
# Notes: None
#
# Copyright: (c) 2005 by the University of Illinois
#
################################################################################

FLUIDHOME=../../..
SRCDIR = $(FLUIDHOME)/utilities/rocflu/init

EXEC= $(FLUIDHOME)/rfluinit

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
vpath PLAG_Mod%.o .:$(FLUIDHOME)/$(BUILDLIB_DIR)
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
main.o: rfluinit.o
UNDERSCORE = $(if $(shell nm rfluinit.o | \
             grep -i rfluinit_),-DFORTRAN_UNDERSCORE,)
LOWERCASE =  $(if $(shell nm rfluinit.o | \
             grep rfluinit),-DFORTRAN_LOWERCASE,)
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
# $Log: Makefile,v $
# Revision 1.4  2008/11/26 15:15:10  mtcampbe
# Port to NCSA Abe cluster for Rocflu.
#
# Revision 1.3  2006/05/05 17:15:28  haselbac
# Added PLAG modules to vpath
#
# Revision 1.2  2005/05/03 03:04:03  haselbac
# Changed to linking with C++
#
# Revision 1.1  2005/04/15 15:05:52  haselbac
# Initial revision
#
################################################################################
