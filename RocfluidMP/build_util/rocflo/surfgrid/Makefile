################################################################################
#
# $Id: Makefile,v 1.6 2003/12/04 03:19:24 haselbac Exp $
#
# Purpose: Makefile for ROCFLO-MP`s surface grid output.
#
# Description: none
#
# Notes: none
#
# Copyright: (c) 2003 by the University of Illinois
#
################################################################################

FLUIDHOME=../../..
SRCDIR = $(FLUIDHOME)/utilities/rocflo/surfgrid

EXEC= $(FLUIDHOME)/rflosurf

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
#   Revision 1.6  2003/12/04 03:19:24  haselbac
#   Modified clean target
#
#   Revision 1.5  2003/11/20 16:40:33  mdbrandy
#   Backing out RocfluidMP changes from 11-17-03
#
#   Revision 1.3  2003/05/15 02:57:00  jblazek
#   Inlined index function.
#
#   Revision 1.2  2003/03/21 23:06:33  haselbac
#   Added LIBINT to deps, made CP robust
#
#   Revision 1.1.1.1  2003/03/20 22:48:04  haselbac
#   Initial revision
#
################################################################################
