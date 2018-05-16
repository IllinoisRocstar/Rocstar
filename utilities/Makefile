SHELL		= /bin/sh
.SUFFIXES:	.o .f90

OSTYPE		= $(shell uname)

#----------------- options for Origin 2000
ifeq "$(OSTYPE)" "IRIX64"
F90		= f90
F90FLAGS        =
endif

#----------------- options for turing cluster
ifeq  "$(OSTYPE)" "Linux"
F90		= f90
F90FLAGS        =
endif

#----------------- options for Blue horizon
ifeq  "$(OSTYPE)" "AIX"
F90		= xlf90
F90FLAGS	= -qsuffix=f=f90
endif

#----------------- options for Blue Pacific
ifeq  "$(OSTYPE)" "PACIFIC"
F90		= xlf90
F90FLAGS	= -qsuffix=f=f90
endif

#----------------- options for Solaris 2.8
ifeq "$(OSTYPE)" "SunOS"
F90             = f90
F90FLAGS        = -xarch=v9
endif

ALL: summary.x select_bbx

summary.x select_bbx:
	$(F90) $(F90FLAGS) $< -o $@

summary.x: summary.f90
select_bbx: select_bbx.f90

clean:
	rm -f *.o summary.x select_bbx
