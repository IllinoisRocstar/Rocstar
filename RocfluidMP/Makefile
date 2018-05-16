################################################################################
#
# $Id: Makefile,v 1.25 2003/05/19 21:18:21 jblazek Exp $
#
# Purpose: Makefile for RocfluidMP
#
# Copyright: (c) 2002 by the University of Illinois
#
################################################################################


### Variables ##################################################################
include Makefile.common


### Targets ####################################################################
all:
	$(MAKE) -C $(BUILDLIB_DIR)
ifndef EXTERN
	$(MAKE) -C $(BUILDSTD_DIR)
	$(MAKE) -C $(BUILDUTIL_DIR)
endif

clean:
	$(MAKE) -C $(BUILDLIB_DIR) $@
	$(MAKE) -C $(BUILDSTD_DIR) $@
	$(MAKE) -C $(BUILDUTIL_DIR) $@

dist:
	$(MAKE) clean
	cd ../
	tar $(TARFLAGS) $(DISTNAME) *
	gzip $(DISTNAME) 

install:
	$(MAKE) -C $(BUILDSTD_DIR) $@
	$(MAKE) -C $(BUILDUTIL_DIR) $@


################ CVS commands
GENX_COMPONENTS = CODING_RULES Makefile Makefile.AIX Makefile.IRIX64 \
        Makefile.Linux Makefile.SunOS Makefile.common README \
        TEMPLATE.F90 genx libflo libfloflu libflu modflo modfloflu \
        modflu rocflo rocflu rocinteract rocpart rocrad rocsmoke rocspecies \
        rocturb rocperi utilities

cvsfwdtag:
	cvs tag -F gen2_5 $(GENX_COMPONENTS)

cvs2genx:
	cvs update -r gen2_5

cvs2main:
	cvs update -A

################################################################################
#
# RCS Revision history:
#
#   $Log: Makefile,v $
#   Revision 1.25  2003/05/19 21:18:21  jblazek
#   Automated switch to 0th-order extrapolation at slip walls and injection boundaries.
#
#   Revision 1.24  2003/04/04 23:17:32  jblazek
#   Added flag for compilation with external driver.
#
#   Revision 1.23  2003/04/04 23:15:53  jblazek
#   Added flag for compilation with external driver.
#
#   Revision 1.22  2003/03/29 03:25:31  wasistho
#   install ROCPERI
#
#   Revision 1.21  2003/03/20 23:46:26  jiao
#   ACH: Merged with Makefile.alone.
#
#   Revision 1.20  2003/03/04 22:12:34  jferry
#   Initial import of Rocinteract
#
#   Revision 1.19  2002/10/23 19:08:57  jiao
#   Removed calcs from gen2
#
#   Revision 1.18  2002/10/15 01:55:54  jblazek
#   Added calcs to cvstag path.
#
#   Revision 1.17  2002/09/27 00:57:07  jblazek
#   Changed makefiles - no makelinks needed.
#
#   Revision 1.16  2002/09/20 22:22:32  jblazek
#   Finalized integration into GenX.
#
################################################################################
