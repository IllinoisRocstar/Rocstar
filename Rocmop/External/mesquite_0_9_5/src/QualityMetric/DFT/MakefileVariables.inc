DFTDIR = $(srcdir)/QualityMetric/DFT

DFTSRC = $(DFTDIR)/I_DFT.cpp \
         $(DFTDIR)/sI_DFT.cpp \
         $(DFTDIR)/RI_DFT.cpp \
         $(DFTDIR)/sRI_DFT.cpp 

DFTHEADERS = \
	$(patsubst $(DFTDIR)/%,$(localincludedir)/%,\
	$(wildcard $(DFTDIR)/*.hpp)) \
	$(patsubst $(DFTDIR)/%,$(localincludedir)/%,\
	$(wildcard $(DFTDIR)/*.hh))

DFTOBJ  = $(patsubst $(DFTDIR)/%.cpp, obj/%.o, $(DFTSRC)) 

ALL_MSQ_OBJ += $(DFTOBJ)
