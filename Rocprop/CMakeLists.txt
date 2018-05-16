cmake_minimum_required (VERSION 2.8)
project(Rocprop)

load_cache(../)
if(NOT ENABLE_MPI)
  add_definitions( -DDUMMY_MPI )
else()
  find_package(MPI REQUIRED)
  add_definitions( -DMPICH_IGNORE_CXX_SEEK )
  include_directories(${MPI_INCLUDE_PATH})
endif()

#INCLUDE(CTest)
find_path(IMPACT_INC com.h HINTS ../IMPACT/COM/include)
find_path(SURF_INC Manifold_2.h HINTS ../IMPACT/SufUtil/include)
find_path(MAP_INC mapbasic.h HINTS ../IMPACT/SurfMap/include)
find_path(SIMPAL_INC Rocblas.h HINTS ../IMPACT/Simpal/include)
set (LIB_SRCS src/Propagation_3.C src/MarkerParticles_3.C src/FaceOffset_3.C src/Rocprop.C
              src/NuLaplacian.C src/cons_diff.C src/AnisotropicSmoothing.C src/quadric_analysis.C
              src/detect_features.C)
set (TEST_SRCS test/proptest.C test/mptest.C test/compare_normals.C test/rotate.C
               test/advectest.C test/detectfeas.C test/PointPropagate.C)
set (ALL_SRCS ${LIB_SRCS} ${TEST_SRCS})

set_source_files_properties(${ALL_SRCS} PROPERTIES COMPILE_FLAGS "-fPIC ${MPI_CXX_COMPILE_FLAGS}" )
#set(TEST_SRCS src/MANTest.C)

# rpath settings
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include_directories(include ${SIMPAL_INC} ${MAP_INC} ${SURF_INC} ${SURF_INC}/../test ${ROCCOM_INC})


if(NOT BUILD_STATIC)
  add_library(Rocprop SHARED ${LIB_SRCS})
else()
  add_library(Rocprop STATIC ${LIB_SRCS})
  add_definitions( -DSTATIC_LINK )
endif()
target_link_libraries(Rocprop ${IMPACT_LIB} ${MPI_CXX_LIBRARIES})

add_executable(proptest test/proptest.C)
add_executable(mptest test/mptest.C)
add_executable(compare_normals test/compare_normals.C)
add_executable(rotate test/rotate.C)
add_executable(advectest test/advectest.C test/PointPropagate.C)
add_executable(detectfeas test/detectfeas.C)
target_link_libraries(proptest Rocprop ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
target_link_libraries(mptest Rocprop ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
target_link_libraries(compare_normals Rocprop ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
target_link_libraries(rotate Rocprop ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
target_link_libraries(advectest Rocprop ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
target_link_libraries(detectfeas Rocprop ${MPI_CXX_LIBRARIES} ${HDF5_LIB})

set_target_properties(Rocprop PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
set_target_properties(proptest PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
set_target_properties(mptest PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
set_target_properties(compare_normals PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
set_target_properties(rotate PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
set_target_properties(advectest PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
set_target_properties(detectfeas PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

install(TARGETS Rocprop RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
