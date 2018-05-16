QUALASSESSDIR = $(srcdir)/QualityAssessor

QUALASSESSSRC = $(QUALASSESSDIR)/QualityAssessor.cpp

QUALASSESSHEADERS = \
	$(patsubst $(QUALASSESSDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QUALASSESSDIR)/*.hpp)) \
	$(patsubst $(QUALASSESSDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QUALASSESSDIR)/*.hh))

QUALASSESSOBJ  = $(patsubst $(QUALASSESSDIR)/%.cpp, obj/%.o, $(QUALASSESSSRC)) 

ALL_MSQ_OBJ += $(QUALASSESSOBJ)
