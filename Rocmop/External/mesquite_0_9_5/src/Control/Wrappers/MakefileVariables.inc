WRAPPERSDIR = $(srcdir)/Control/Wrappers

WRAPPERSSRC = $(WRAPPERSDIR)/ShapeImprovementWrapper.cpp

WRAPPERSHEADERS = \
	$(patsubst $(WRAPPERSDIR)/%,$(localincludedir)/%,\
	$(wildcard $(WRAPPERSDIR)/*.hpp)) \
	$(patsubst $(WRAPPERSDIR)/%,$(localincludedir)/%,\
	$(wildcard $(WRAPPERSDIR)/*.hh))

WRAPPERSOBJ  = $(patsubst $(WRAPPERSDIR)/%.cpp, obj/%.o, $(WRAPPERSSRC))

ALL_MSQ_OBJ += $(WRAPPERSOBJ)
