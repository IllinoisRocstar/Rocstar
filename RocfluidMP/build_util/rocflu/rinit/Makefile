################################################################################
#
# $Id: Makefile,v 1.1 2004/03/19 21:55:16 haselbac Exp $
#
# Purpose: Makefile for ROCFLU-MP's reinitializor.
#
# Description: none
#
# Notes: none
#
# Copyright: (c) 2004 by the University of Illinois
#
################################################################################

FLUIDHOME=../../..
SRCDIR = $(FLUIDHOME)/utilities/rocflu/rinit

EXEC= $(FLUIDHOME)/rflurinit

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
# $Log: Makefile,v $
# Revision 1.1  2004/03/19 21:55:16  haselbac
# Initial revision
#
################################################################################
