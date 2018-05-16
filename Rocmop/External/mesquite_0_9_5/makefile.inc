# locallibdir - where the compiled libraries will be linked
locallibdir = ./lib
# locallibdir - where all the include files will be (locally) copied.
localincludedir = ./includeLinks
# localobjdir - where all the object files will be written
localobjdir = ./obj
# file where to write dependencies - will be included in this Makefile
dependenciesfile = make.dependencies
# template_dir - where template instantiations are stored.
#                Only applies to some platforms.
template_dir = 

#include Makefile.customize

OUTPUT_OPTION = -o $@ # Use the -o option when generating .o files.
INCLUDE = ${SYSTEM_INCLUDE} -I$(SRCDIR)/include -I$(localincludedir) -I./TSTT-interface

# Add to this in each subdirectory.  It's the list of files
# that get added to the mesquite library.
ALL_MSQ_OBJ = 

# Make sure that "default" is the first target
# default: depend all

# List all desired module names explicitly below
srcdir = $(SRCDIR)/src
MODULENAMES :=  Mesh \
		Control \
		Control/Wrappers \
		Misc \
                ObjectiveFunction \
		QualityAssessor \
		QualityImprover \
		QualityImprover/TopologyModifier \
		QualityImprover/VertexMover \
		QualityImprover/VertexMover/ConjugateGradient \
		QualityImprover/VertexMover/LaplacianSmoothers \
		QualityImprover/VertexMover/NonSmoothSteepestDescent \
		QualityImprover/VertexMover/SteepestDescent \
		QualityImprover/VertexMover/FeasibleNewton \
		QualityImprover/VertexMover/Randomize \
		QualityMetric \
		QualityMetric/Shape \
		QualityMetric/Smoothness \
		QualityMetric/Untangle \
		QualityMetric/Volume \
		QualityMetric/DFT \
                TargetCalculator \
		lib 

testdir = testSuite
TESTNAMES := test_1\
             laplacian_test\
             laplacian_wrapper_test\
             untangle_test\
             analytical_grad_3D\
             algorithm_test\
             ActiveSetTest\
             simple_hybrid_test\
             wrapper_tests\
             tutorial\
             test_DFT \
             escobar \
             Guides843 
###             random_test\
###             convert\
###             shape_wrapper_test

# ************ inclusion of all the modules 
# ************ MakefileVariables.inc  and
# ************ MakefileTargets.inc

# Add location to the beginning of each (i.e. './')
MODULES := $(wildcard $(patsubst %, $(srcdir)/%, $(MODULENAMES)))

# Generate a list of module makefiles actually present. 
MODULEMAKEFILES := $(wildcard $(patsubst %, %/MakefileVariables.inc,\
                    $(MODULES)))

# the module directories are added to the include path
# not in use - using links setup includeLinks instead 
#INCLUDES += $(patsubst %, -I%, $(MODULES))

# add all .cpp and .cc files to the list of sources. This
# list will be sent to makedepend to automatically generate
# dependancies.
ALLSRC := $(foreach MODULE, $(MODULES),\
	 $(wildcard $(MODULE)/*.cpp $(srcdir)/$(MODULE)/*.cc)) 

# now include the module makefiles (if there are any)
ifdef MODULEMAKEFILES
include $(MODULEMAKEFILES)
endif

# include any module targets, such as executables, that cannot
# be included in the module makefiles because they require the
# value of variables from other modules (make expands variables
# in the dependancy lists as soon as the target is read)
#MODULETARGETFILES = $(wildcard $(srcdir)/*/MakefileTargets.inc)

MODULETARGETFILES := $(wildcard $(patsubst %, \
	%/MakefileTargets.inc, $(MODULES)))

ifdef MODULETARGETFILES
include $(MODULETARGETFILES)
endif
# ************

