FEASNEWTDIR = $(srcdir)/QualityImprover/VertexMover/FeasibleNewton

FEASNEWTSRC = $(FEASNEWTDIR)/FeasibleNewton.cpp

FEASNEWTHEADERS = \
	$(patsubst $(FEASNEWTDIR)/%,$(localincludedir)/%,\
	$(wildcard $(FEASNEWTDIR)/*.hpp)) \
	$(patsubst $(FEASNEWTDIR)/%,$(localincludedir)/%,\
	$(wildcard $(FEASNEWTDIR)/*.hh))

FEASNEWTOBJ  = $(patsubst $(FEASNEWTDIR)/%.cpp, obj/%.o, $(FEASNEWTSRC)) 

ALL_MSQ_OBJ += $(FEASNEWTOBJ)
