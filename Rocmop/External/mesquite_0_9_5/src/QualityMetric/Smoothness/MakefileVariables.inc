QMSMOOTHNESSDIR = $(srcdir)/QualityMetric/Smoothness

QMSMOOTHNESSSRC = $(QMSMOOTHNESSDIR)/EdgeLengthQualityMetric.cpp \
                  $(QMSMOOTHNESSDIR)/EdgeLengthRangeQualityMetric.cpp \
                  $(QMSMOOTHNESSDIR)/ASMQualityMetric.cpp \


QMSMOOTHNESSHEADERS = \
	$(patsubst $(QMSMOOTHNESSDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMSMOOTHNESSDIR)/*.hpp)) \
	$(patsubst $(QMSMOOTHNESSDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMSMOOTHNESSDIR)/*.hh))

QMSMOOTHNESSOBJ  = $(patsubst $(QMSMOOTHNESSDIR)/%.cpp, obj/%.o, $(QMSMOOTHNESSSRC)) 

ALL_MSQ_OBJ += $(QMSMOOTHNESSOBJ)
