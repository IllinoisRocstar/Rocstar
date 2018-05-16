MISCDIR = $(srcdir)/Misc

MISCSRCFILES = Exponent.cpp \
               FileTokenizer.cpp \
               MeshTransform.cpp \
               MesquiteVersion.cpp \
               MsqDebug.cpp \
               MsqError.cpp \
               MsqFPE.cpp \
               MsqHessian.cpp \
               MsqInterrupt.cpp \
               MsqTimer.cpp \
               PlanarDomain.cpp \
               SphericalDomain.cpp \
               Vector3D.cpp

MISCSRC = $(MISCSRCFILES:%=$(MISCDIR)/%)

MISCHEADERS = \
	$(patsubst $(MISCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(MISCDIR)/*.hpp)) \
	$(patsubst $(MISCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(MISCDIR)/*.hh))

MISCOBJ  = $(patsubst $(MISCDIR)/%.cpp, obj/%.o, $(MISCSRC))

ALL_MSQ_OBJ += $(MISCOBJ)
