NONSMOOTHSTEEPDESCDIR = $(srcdir)/QualityImprover/VertexMover/NonSmoothSteepestDescent

NONSMOOTHSTEEPDESCSRC = $(NONSMOOTHSTEEPDESCDIR)/NonSmoothSteepestDescent.cpp

NONSMOOTHSTEEPDESCHEADERS = \
	$(patsubst $(NONSMOOTHSTEEPDESCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(NONSMOOTHSTEEPDESCDIR)/*.hpp)) \
	$(patsubst $(NONSMOOTHSTEEPDESCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(NONSMOOTHSTEEPDESCDIR)/*.hh))

NONSMOOTHSTEEPDESCOBJ  = $(patsubst $(NONSMOOTHSTEEPDESCDIR)/%.cpp, obj/%.o,$(NONSMOOTHSTEEPDESCSRC)) 

ALL_MSQ_OBJ += $(NONSMOOTHSTEEPDESCOBJ)
