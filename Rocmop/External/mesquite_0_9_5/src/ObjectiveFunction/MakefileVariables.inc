OBJFUNCDIR = $(srcdir)/ObjectiveFunction

OBJFUNCSRC = $(OBJFUNCDIR)/CompositeOFAdd.cpp\
	       $(OBJFUNCDIR)/CompositeOFMultiply.cpp \
	       $(OBJFUNCDIR)/CompositeOFScalarAdd.cpp \
	       $(OBJFUNCDIR)/CompositeOFScalarMultiply.cpp \
	       $(OBJFUNCDIR)/LInfTemplate.cpp \
	       $(OBJFUNCDIR)/LPtoPTemplate.cpp \
	       $(OBJFUNCDIR)/MaxTemplate.cpp \
	       $(OBJFUNCDIR)/ObjectiveFunction.cpp 

OBJFUNCHEADERS = \
	$(patsubst $(OBJFUNCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(OBJFUNCDIR)/*.hpp)) \
	$(patsubst $(OBJFUNCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(OBJFUNCDIR)/*.hh))

OBJFUNCOBJ  = $(patsubst $(OBJFUNCDIR)/%.cpp, obj/%.o, $(OBJFUNCSRC))

ALL_MSQ_OBJ += $(OBJFUNCOBJ)
