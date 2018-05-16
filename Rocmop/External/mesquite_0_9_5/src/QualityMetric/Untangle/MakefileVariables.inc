QMUNTANGLEDIR = $(srcdir)/QualityMetric/Untangle

QMUNTANGLESRC = $(QMUNTANGLEDIR)/UntangleBetaQualityMetric.cpp

QMUNTANGLEHEADERS = \
	$(patsubst $(QMUNTANGLEDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMUNTANGLEDIR)/*.hpp)) \
	$(patsubst $(QMUNTANGLEDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMUNTANGLEDIR)/*.hh))

QMUNTANGLEOBJ  = $(patsubst $(QMUNTANGLEDIR)/%.cpp, obj/%.o, $(QMUNTANGLESRC)) 

ALL_MSQ_OBJ += $(QMUNTANGLEOBJ)
