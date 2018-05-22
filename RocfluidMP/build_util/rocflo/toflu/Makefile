################################################################################
#
# $Id: Makefile,v 1.1.1.1 2004/08/17 01:46:05 wasistho Exp $
#
# Purpose: Makefile for FLO2FLU conversion.
#
# Description: none
#
# Notes: none
#
# Copyright: (c) 2003 by the University of Illinois
#
################################################################################

FLUIDHOME=../../..
SRCDIR = $(FLUIDHOME)/utilities/rocflo/toflu

EXEC= $(FLUIDHOME)/rflo2flu

LIBINT=$(FLUIDHOME)/$(BUILDLIB_DIR)/libflo.a

MODDIR=$(FLUIDHOME)/$(BUILDLIB_DIR)
SRCF77=
SRCF90=
SYSTEM=$(shell uname -s)

FPPFLAGS+= -I$(FLUIDHOME)/modflo

### Include common and dependent makefiles #####################################
include $(FLUIDHOME)/Makefile.common
include $(FLUIDHOME)/Makefile.$(SYSTEM)
LIBEXT= $(LIBSYS)

### Define SRCF90 and SRCF77 ###################################################
include $(SRCDIR)/$(FILELIST)

vpath %.F90 $(SRCDIR)
vpath Mod%.o .:$(FLUIDHOME)/$(BUILDLIB_DIR)
vpath RFLO_Mod%.o .:$(FLUIDHOME)/$(BUILDLIB_DIR)

include $(FLUIDHOME)/Makefile.dep

### Rules ######################################################################
OBJF77=$(SRCF77:.f=.o)
OBJF90=$(SRCF90:.F90=.o)


### Targets ####################################################################
.PHONY:	all clean

all:    $(EXEC)

$(EXEC): $(OBJF77) $(OBJF90) $(LIBINT)                     
	$(COMP) $(OBJF77) $(OBJF90) -o $(EXEC) $(LDFLAGS) $(LIBINT) $(LIBEXT)

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
#   Revision 1.1.1.1  2004/08/17 01:46:05  wasistho
#   initial import
#
#
################################################################################
