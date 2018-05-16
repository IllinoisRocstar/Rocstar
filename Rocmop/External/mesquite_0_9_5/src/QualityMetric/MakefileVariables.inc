QUALMETRICDIR = $(srcdir)/QualityMetric

QUALMETRICSRC = $(QUALMETRICDIR)/AddQualityMetric.cpp \
		$(QUALMETRICDIR)/MultiplyQualityMetric.cpp \
                $(QUALMETRICDIR)/PowerQualityMetric.cpp \
                $(QUALMETRICDIR)/ScalarAddQualityMetric.cpp \
                $(QUALMETRICDIR)/ScalarMultiplyQualityMetric.cpp \
                $(QUALMETRICDIR)/QualityMetric.cpp

QUALMETRICHEADERS = \
	$(patsubst $(QUALMETRICDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QUALMETRICDIR)/*.hpp)) \
	$(patsubst $(QUALMETRICDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QUALMETRICDIR)/*.hh))

QUALMETRICOBJ  = $(patsubst $(QUALMETRICDIR)/%.cpp, obj/%.o, $(QUALMETRICSRC)) 

ALL_MSQ_OBJ += $(QUALMETRICOBJ)
