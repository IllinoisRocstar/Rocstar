################################################################################
#
# $Id: Makefile,v 1.14 2008/11/26 15:15:23 mtcampbe Exp $
#
# Purpose: Makefile for ROCFLU-MP's postprocessor.
#
# Description: None.
#
# Notes: None.
#
# Copyright: (c) 2001-2005 by the University of Illinois
#
################################################################################

FLUIDHOME=../../..
SRCDIR = $(FLUIDHOME)/utilities/rocflu/post
EXEC= $(FLUIDHOME)/rflupost

LIBINT=$(FLUIDHOME)/$(BUILDLIB_DIR)/libflu.a

MODDIR=$(FLUIDHOME)/$(BUILDLIB_DIR)
SRCF77=
SRCF90=
SYSTEM=$(shell uname -s)
RFLUPOST=1

### Include common and dependent makefiles #####################################
include $(FLUIDHOME)/Makefile.common
include $(FLUIDHOME)/Makefile.$(SYSTEM)
LIBEXT= $(LIBTEC) $(LIBSYS)

### Define SRCF90 and SRCF77 ###################################################
include $(SRCDIR)/$(FILELIST)

vpath %.F90 $(SRCDIR)
vpath Mod%.o        .:$(FLUIDHOME)/$(BUILDLIB_DIR)
vpath RFLU_Mod%.o   .:$(FLUIDHOME)/$(BUILDLIB_DIR)
vpath %.C $(SRCDIR)

ifdef PLAG
  vpath PLAG_Mod%.o  $(FLUIDHOME)/$(BUILDLIB_DIR)
endif

ifdef TURB
  vpath TURB_Mod%.o  $(FLUIDHOME)/$(BUILDLIB_DIR)
endif

include $(FLUIDHOME)/Makefile.dep

### Rules ######################################################################
OBJF77=$(SRCF77:.f=.o)
OBJF90=$(SRCF90:.F90=.o)
OBJCXX=$(SRCCXX:.C=.o)

### Targets ####################################################################
.PHONY:	all clean

all:    $(EXEC)

# Automatically determine the Fortran name mangling
main.o: rflupost.o
UNDERSCORE = $(if $(shell nm rflupost.o | \
             grep -i rflupost_),-DFORTRAN_UNDERSCORE,)
LOWERCASE =  $(if $(shell nm rflupost.o | \
             grep rflupost),-DFORTRAN_LOWERCASE,)
CXXFLAGS  += $(LOWERCASE) $(UNDERSCORE)

$(EXEC): $(OBJF77) $(OBJF90) $(LIBINT) $(OBJCXX)                    
	$(LINKER) $(OBJF77) $(OBJF90) $(OBJCXX) -o $(EXEC) $(LDFLAGS) $(LIBINT) $(LIBEXT) $(CXXLINKFLAGS)

clean:  
	-$(RM) $(EXEC)
	-$(RM) $(MODULES)
	-$(RM) $(RIFFRAFF)

install:
	-$(CP) $(EXEC) $(BINDIR)/.


################################################################################
#
# RCS Revision history:
#
# $Log: Makefile,v $
# Revision 1.14  2008/11/26 15:15:23  mtcampbe
# Port to NCSA Abe cluster for Rocflu.
#
# Revision 1.13  2006/01/12 09:45:36  wasistho
# added access to TURB_Mod
#
# Revision 1.12  2005/07/05 16:15:48  haselbac
# Removed CHARM ifdef
#
# Revision 1.11  2005/05/05 18:31:36  haselbac
# Renamed TECINT variable to RFLUPOST
#
# Revision 1.10  2005/05/03 20:37:08  haselbac
# Added TECINT to work around TECPLOT linking conflict problem, cosmetics
#                                   
# Revision 1.9  2005/05/03 03:04:04  haselbac          
# Changed to linking with C++                          
#
# Revision 1.8  2005/03/24 18:41:12  haselbac          
# Reversed order of libraries - system should go last  
#
# Revision 1.7  2004/03/05 23:22:24  haselbac          
# Added vpath for PLAG                                 
#
# Revision 1.6  2003/12/04 03:19:30  haselbac          
# Modified clean target                                
#
# Revision 1.5  2003/04/11 15:15:18  jiao              
# Changed RFLU_CHARM back to CHARM.                    
#
# Revision 1.4  2003/03/25 19:10:05  haselbac          
# Added - to RM                                        
#
# Revision 1.3  2003/03/22 00:29:29  haselbac          
# Removed old footer entries                           
#
# Revision 1.2  2003/03/21 23:07:17  haselbac          
# Added LIBINT to deps, made CP robust                 
#
# Revision 1.1.1.1  2003/03/20 22:48:04  haselbac      
# Initial revision                                     
#
################################################################################
