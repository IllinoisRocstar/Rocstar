################################################################################
#
# $Id: Makefile,v 1.12 2005/05/03 03:03:05 haselbac Exp $
#
# Purpose: Makefile for RocfluidMP libraries
#
# Copyright: (c) 2003 by the University of Illinois
#
################################################################################

ifdef RFLO
  LIBFLUID= libflo.a
endif
ifdef RFLU
  LIBFLUID= libflu.a
endif

FPPFLAGS=
SRCF77=
SRCF90=
SYSTEM= $(shell uname -s)
FLUIDHOME= ../
FPPFLAGS+= -I../modflo

### Include common and dependent makefiles #####################################
include ../Makefile.common
include ../Makefile.$(SYSTEM)

### vpath definitions ##########################################################
vpath Mod%.F90 		../modfloflu:../modflo

ifdef RFLO
  vpath RFLO_Mod%.F90	../modflo
  vpath RFLO_%.F90	../libflo:../rocflo
  vpath %.F90		../libfloflu:../rocflo
endif
ifdef RFLU
  vpath RFLU_Mod%.F90	../modflu
  vpath RFLU_%.F90	../libflu
  vpath %.F90           ../libfloflu:../rocflu
  vpath %.f		../libflu
  vpath %.C		../libflu
endif

ifdef INRT
  vpath INRT_Mod%.F90   ../rocinteract
  vpath INRT%.F90       ../rocinteract
endif

ifdef PEUL
  vpath PEUL_Mod%.F90   ../rocsmoke
  vpath PEUL%.F90       ../rocsmoke
endif

ifdef PLAG
  vpath PLAG_Mod%.F90   ../rocpart
  vpath PLAG%.F90       ../rocpart
endif

ifdef RADI
  vpath RADI_Mod%.F90   ../rocrad
  vpath RADI%.F90       ../rocrad
endif

ifdef SPEC
  vpath SPEC_Mod%.F90        ../rocspecies
  vpath SPEC_RFLU_Mod%.F90   ../rocspecies
  vpath SPEC%.F90            ../rocspecies
endif

ifdef TURB
  vpath TURB_Mod%.F90   ../rocturb
  vpath TURB%.F90       ../rocturb
endif

ifdef PERI
  vpath PERI_Mod%.F90   ../rocperi
  vpath PERI%.F90       ../rocperi
endif

### Define SRCF90 and SRCF77 ###################################################
include ../modfloflu/$(FILELIST)

ifdef RFLO
  include ../modflo/$(FILELIST)
endif
ifdef RFLU
  include ../modflu/$(FILELIST)
endif

include ../libfloflu/$(FILELIST)

ifdef RFLO
  include ../libflo/$(FILELIST)
  include ../rocflo/$(FILELIST)
endif
ifdef RFLU
  include ../libflu/$(FILELIST)
  include ../rocflu/$(FILELIST)
endif

ifdef INRT
  include ../rocinteract/$(FILELIST)
endif

ifdef PEUL
  include ../rocsmoke/$(FILELIST)
endif

ifdef PLAG
  include ../rocpart/$(FILELIST)
endif

ifdef RADI
  include ../rocrad/$(FILELIST)
endif

ifdef SPEC
  include ../rocspecies/$(FILELIST)
endif

ifdef TURB
  include ../rocturb/$(FILELIST)
endif

ifdef PERI
  include ../rocperi/$(FILELIST)
endif

### Rules ######################################################################
OBJF77= $(SRCF77:.f=.o)
OBJF90= $(SRCF90:.F90=.o)
OBJCXX= $(SRCCXX:.C=.o)

### Targets ####################################################################
all:    $(LIBFLUID)

$(LIBFLUID): $(OBJF77) $(OBJF90) $(OBJCXX)
	$(AR) $(ARFLAGS) $@ $?                        

clean:  
	$(FIND) "$(RIFFRAFF)" | $(XARGS) $(RM)
	$(RM) $(MODULES)
	$(RM) $(LIBFLUID)

include ../Makefile.dep

################################################################################
#
# RCS Revision history:
#
#   $Log: Makefile,v $
#   Revision 1.12  2005/05/03 03:03:05  haselbac
#   Added compilation of C++ files
#
#   Revision 1.11  2004/04/01 21:23:16  haselbac
#   Added vpath for SPEC RFLU modules
#
#   Revision 1.10  2003/12/04 03:19:11  haselbac
#   Modified clean target
#
#   Revision 1.9  2003/11/20 16:40:32  mdbrandy
#   Backing out RocfluidMP changes from 11-17-03
#
#   Revision 1.5  2003/07/03 21:48:44  jblazek
#   Implemented dual-time stepping.
#
#   Revision 1.4  2003/05/15 02:56:59  jblazek
#   Inlined index function.
#
#   Revision 1.3  2003/03/29 03:25:57  wasistho
#   install ROCPERI
#
#   Revision 1.2  2003/03/21 17:09:44  haselbac
#   Trimmed footer
#
#   Revision 1.1.1.1  2003/03/20 22:45:06  haselbac
#   Initial revision
#
################################################################################
