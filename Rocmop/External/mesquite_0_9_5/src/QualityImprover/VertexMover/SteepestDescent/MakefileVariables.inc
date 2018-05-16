STEEPDESCDIR = $(srcdir)/QualityImprover/VertexMover/SteepestDescent

STEEPDESCSRC = $(STEEPDESCDIR)/SteepestDescent.cpp

STEEPDESCHEADERS = \
	$(patsubst $(STEEPDESCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(STEEPDESCDIR)/*.hpp)) \
	$(patsubst $(STEEPDESCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(STEEPDESCDIR)/*.hh))

STEEPDESCOBJ  = $(patsubst $(STEEPDESCDIR)/%.cpp, obj/%.o, $(STEEPDESCSRC)) 

ALL_MSQ_OBJ += $(STEEPDESCOBJ)
