TOPOMODIFDIR = $(srcdir)/QualityImprover/TopologyModifier

TOPOMODIFSRC = $(TOPOMODIFDIR)/TopologyModifier.cpp

TOPOMODIFHEADERS = \
	$(patsubst $(TOPOMODIFDIR)/%,$(localincludedir)/%,\
	$(wildcard $(TOPOMODIFDIR)/*.hpp)) \
	$(patsubst $(TOPOMODIFDIR)/%,$(localincludedir)/%,\
	$(wildcard $(TOPOMODIFDIR)/*.hh))

TOPOMODIFOBJ  = $(patsubst $(TOPOMODIFDIR)/%.cpp, obj/%.o, $(TOPOMODIFSRC)) 

ALL_MSQ_OBJ += $(TOPOMODIFOBJ)
