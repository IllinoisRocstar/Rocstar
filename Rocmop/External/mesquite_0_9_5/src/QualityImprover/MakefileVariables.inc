QUALIMPROVERDIR = $(srcdir)/QualityImprover

QUALIMPROVERSRC = $(QUALIMPROVERDIR)/QualityImprover.cpp 

QUALIMPROVERHEADERS = \
	$(patsubst $(QUALIMPROVERDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QUALIMPROVERDIR)/*.hpp)) \
	$(patsubst $(QUALIMPROVERDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QUALIMPROVERDIR)/*.hh))

QUALIMPROVEROBJ  = $(patsubst $(QUALIMPROVERDIR)/%.cpp, obj/%.o, $(QUALIMPROVERSRC))

ALL_MSQ_OBJ += $(QUALIMPROVEROBJ)
