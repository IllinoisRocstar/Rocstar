RANDSMOOTHDIR = $(srcdir)/QualityImprover/VertexMover/Randomize

RANDSMOOTHSRC = $(RANDSMOOTHDIR)/Randomize.cpp 

RANDSMOOTHHEADERS = \
	$(patsubst $(RANDSMOOTHDIR)/%,$(localincludedir)/%,\
	$(wildcard $(RANDSMOOTHDIR)/*.hpp)) \
	$(patsubst $(RANDSMOOTHDIR)/%,$(localincludedir)/%,\
	$(wildcard $(RANDSMOOTHDIR)/*.hh))

RANDSMOOTHOBJ  = $(patsubst $(RANDSMOOTHDIR)/%.cpp, obj/%.o, $(RANDSMOOTHSRC)) 

ALL_MSQ_OBJ += $(RANDSMOOTHOBJ)
