CONJGRADDIR = $(srcdir)/QualityImprover/VertexMover/ConjugateGradient

CONJGRADSRC = $(CONJGRADDIR)/ConjugateGradient.cpp

CONJGRADHEADERS = \
	$(patsubst $(CONJGRADDIR)/%,$(localincludedir)/%,\
	$(wildcard $(CONJGRADDIR)/*.hpp)) \
	$(patsubst $(CONJGRADDIR)/%,$(localincludedir)/%,\
	$(wildcard $(CONJGRADDIR)/*.hh))

CONJGRADOBJ  = $(patsubst $(CONJGRADDIR)/%.cpp, obj/%.o, $(CONJGRADSRC))

ALL_MSQ_OBJ += $(CONJGRADOBJ)
