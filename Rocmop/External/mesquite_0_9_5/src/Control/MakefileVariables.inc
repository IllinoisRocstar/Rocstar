CONTROLDIR = $(srcdir)/Control

CONTROLSRC = $(CONTROLDIR)/InstructionQueue.cpp \
	     $(CONTROLDIR)/TerminationCriterion.cpp 

CONTROLHEADERS = \
	$(patsubst $(CONTROLDIR)/%,$(localincludedir)/%,\
	$(wildcard $(CONTROLDIR)/*.hpp)) \
	$(patsubst $(CONTROLDIR)/%,$(localincludedir)/%,\
	$(wildcard $(CONTROLDIR)/*.hh))

CONTROLOBJ  = $(patsubst $(CONTROLDIR)/%.cpp, obj/%.o, $(CONTROLSRC))

ALL_MSQ_OBJ += $(CONTROLOBJ)
