################################################################################
#
# $Id: Makefile,v 1.3 2008/11/26 15:15:16 mtcampbe Exp $
#
# Purpose: Makefile for ROCFLU-MP's partitioning utility.
#
# Description: None.
#
# Notes: None.
#
# Copyright: (c) 2005 by the University of Illinois
#
################################################################################

FLUIDHOME=../../..
SRCDIR = $(FLUIDHOME)/utilities/rocflu/part

EXEC= $(FLUIDHOME)/rflupart

LIBINT=$(FLUIDHOME)/$(BUILDLIB_DIR)/libflu.a

MODDIR=$(FLUIDHOME)/$(BUILDLIB_DIR)
SRCF77=
SRCF90=
SYSTEM=$(shell uname -s)

### Include common and dependent makefiles #####################################
include $(FLUIDHOME)/Makefile.common
include $(FLUIDHOME)/Makefile.$(SYSTEM)
LIBEXT= $(LIBSYS) $(LIBMET)

### Define SRCF90 and SRCF77 ###################################################
include $(SRCDIR)/$(FILELIST)

vpath %.C $(SRCDIR)
vpath %.F90 $(SRCDIR)
vpath Mod%.o .:$(FLUIDHOME)/$(BUILDLIB_DIR)
vpath RFLU_Mod%.o .:$(FLUIDHOME)/$(BUILDLIB_DIR)

include $(FLUIDHOME)/Makefile.dep

### Rules ######################################################################
OBJF77=$(SRCF77:.f=.o)
OBJF90=$(SRCF90:.F90=.o)
OBJCXX=$(SRCCXX:.C=.o)

### Targets ####################################################################
.PHONY:	all clean

all:    $(EXEC)

# Automatically determine the Fortran name mangling
main.o: rflupart.o
UNDERSCORE = $(if $(shell nm rflupart.o | \
             grep -i rflupart_),-DFORTRAN_UNDERSCORE,)
LOWERCASE =  $(if $(shell nm rflupart.o | \
             grep rflupart),-DFORTRAN_LOWERCASE,)
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
# Revision 1.3  2008/11/26 15:15:16  mtcampbe
# Port to NCSA Abe cluster for Rocflu.
#
# Revision 1.2  2005/05/03 03:04:03  haselbac
# Changed to linking with C++
#
# Revision 1.1  2005/04/15 15:05:56  haselbac
# Initial revision
# 
################################################################################
