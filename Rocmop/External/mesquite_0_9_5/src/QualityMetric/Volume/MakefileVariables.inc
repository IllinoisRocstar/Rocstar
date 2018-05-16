QMVOLUMEDIR = $(srcdir)/QualityMetric/Volume

QMVOLUMESRC = $(QMVOLUMEDIR)/CornerJacobianQualityMetric.cpp\
              $(QMVOLUMEDIR)/LocalSizeQualityMetric.cpp\

QMVOLUMEHEADERS = \
	$(patsubst $(QMVOLUMEDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMVOLUMEDIR)/*.hpp)) \
	$(patsubst $(QMVOLUMEDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMVOLUMEDIR)/*.hh))

QMVOLUMEOBJ  = $(patsubst $(QMVOLUMEDIR)/%.cpp, obj/%.o, $(QMVOLUMESRC))

ALL_MSQ_OBJ += $(QMVOLUMEOBJ)
