TARGCALCDIR = $(srcdir)/TargetCalculator

TARGCALCSRC = $(TARGCALCDIR)/TargetCalculator.cpp \
   $(TARGCALCDIR)/LVQDTargetCalculator.cpp \
   $(TARGCALCDIR)/WTargetCalculator.cpp 

TARGCALCHEADERS = \
	$(patsubst $(TARGCALCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(TARGCALCDIR)/*.hpp)) \
	$(patsubst $(TARGCALCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(TARGCALCDIR)/*.hh))

TARGCALCOBJ  = $(patsubst $(TARGCALCDIR)/%.cpp, obj/%.o, $(TARGCALCSRC))

ALL_MSQ_OBJ += $(TARGCALCOBJ)
