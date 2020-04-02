cmake_minimum_required(VERSION 3.0)

set(ROCMOP_SRCS
    src/Rocmop_1.C
    src/smooth_mesquite_1.C
    src/smooth_mesquite_ng_1.C
    src/smooth_medial_1.C)
set(ROCMOP2_SRCS
    src/Rocmop_2.C
    src/smooth_mesquite_2.C
    src/smooth_mesquite_ng_2.C
    src/smooth_medial_2.C)

set(LIB_SRCS
    src/geometry.C
    src/Matrix.C
    src/Algebraic_Metrics_2.C
    src/Algebraic_Metrics_3.C
    src/Geometric_Metrics_2.C
    src/Geometric_Metrics_3.C
    src/PN_patch.C)

if(ROCMOP2)
  set(LIB_SRCS ${LIB_SRCS} ${ROCMOP2_SRCS} src/MesqPane_95.C)
else()
  set(LIB_SRCS ${LIB_SRCS} ${ROCMOP_SRCS} src/MesqPane_95.C)
endif()

set(UTIL_SRCS util/chkPconnGRecv.cc util/RocMeshData.cc)

add_subdirectory(External/mesquite_0_9_5)
add_definitions(-DMESQUITE)

if(ENABLE_DEBUG)
    add_definitions(-DMSQ_ENABLE_DEBUG)
endif(ENABLE_DEBUG)

add_library(Rocmop ${LIB_SRCS})
target_include_directories(Rocmop
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/rocstar>)
target_link_libraries(Rocmop IMPACT::SurfUtil mesquite_95)
set_target_properties(Rocmop PROPERTIES VERSION ${ROCSTAR_VERSION}
        SOVERSION ${ROCSTAR_MAJOR_VERSION})

add_executable(chkPconnGRecv util/chkPconnGRecv.cc util/RocMeshData.cc)
add_executable(add_aspect_ratios test/add_aspect_ratios.C)
add_executable(build_meshes test/build_meshes.C)
add_executable(metric_demo test/metric_demo.C)
add_executable(smooth_volume test/smooth_volume.C)
target_link_libraries(chkPconnGRecv Rocmop)
target_link_libraries(add_aspect_ratios Rocmop)
target_link_libraries(build_meshes Rocmop)
target_link_libraries(metric_demo Rocmop)
target_link_libraries(smooth_volume Rocmop)

install(TARGETS Rocmop
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
