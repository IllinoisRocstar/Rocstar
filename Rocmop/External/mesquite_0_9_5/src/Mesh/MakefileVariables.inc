MESHDIR = $(srcdir)/Mesh

MESHSRCFILES = CornerTag.cpp \
               MeshSet.cpp \
               MsqMeshEntity.cpp \
               MsqVertex.cpp \
               PatchData.cpp \
               MeshImpl.cpp \
               MeshImplData.cpp \
               MeshImplTags.cpp \
               TopologyInfo.cpp

ifeq ($(MSQ_TSTT_MESH),yes)
  MESHSRCFILES += MeshTSTT.cpp
endif

ifeq ($(MSQ_TSTT_GEOM),yes)
  MESHSRCFILES += GeomTSTT.cpp
endif

MESHSRC = $(MESHSRCFILES:%=$(MESHDIR)/%)

MESHHEADERS = \
	$(patsubst $(MESHDIR)/%,$(localincludedir)/%,\
	$(wildcard $(MESHDIR)/*.hpp)) \
	$(patsubst $(MESHDIR)/%,$(localincludedir)/%,\
	$(wildcard $(MESHDIR)/*.hh))

MESHOBJ  = $(patsubst $(MESHDIR)/%.cpp, obj/%.o,$(MESHSRC))

ALL_MSQ_OBJ += $(MESHOBJ)
