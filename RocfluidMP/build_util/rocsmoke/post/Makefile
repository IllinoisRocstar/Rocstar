################################################################################
#
# $Id: Makefile,v 1.2 2004/03/25 21:13:25 jferry Exp $
#
# Purpose: Makefile for ROCFLO-MP's postprocessor.
#
# Description: none
#
# Notes: none
#
# Copyright: (c) 2003 by the University of Illinois
#
################################################################################

FLUIDHOME=../../..
SRCDIR = $(FLUIDHOME)/utilities/rocsmoke/post

EXEC= $(FLUIDHOME)/peulpost

LIBINT=$(FLUIDHOME)/$(BUILDLIB_DIR)/libflo.a

MODDIR=$(FLUIDHOME)/$(BUILDLIB_DIR)
SRCF77=
SRCF90=
SYSTEM=$(shell uname -s)

FPPFLAGS+= -I$(FLUIDHOME)/modflo

### Include common and dependent makefiles #####################################
include $(FLUIDHOME)/Makefile.common
include $(FLUIDHOME)/Makefile.$(SYSTEM)
LIBEXT= $(LIBSYS) $(LIBTEC)

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
#   Revision 1.2  2004/03/25 21:13:25  jferry
#   added removal of module files to cleanup
#
#   Revision 1.1  2003/09/25 15:40:22  jferry
#   Implented Rocsmoke post-processing
#
################################################################################
