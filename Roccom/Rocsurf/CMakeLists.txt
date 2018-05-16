cmake_minimum_required (VERSION 2.8)
project(Rocsurf)
#ENABLE_LANGUAGE( CXX Fortran )
#include(FortranCInterface)
#FortranCInterface_HEADER(${PROJECT_SOURCE_DIR}/include/FC.h MACRO_NAMESPACE "FC_")
#FortranCInterface_VERIFY(CXX QUIET)
#set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single directory for all executables.")
#set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Single directory for all libraries and archives.")
#mark_as_advanced (LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)
load_cache(../)
IF(NOT ENABLE_MPI)
  add_definitions( -DDUMMY_MPI )
ELSE()
  FIND_PACKAGE(MPI REQUIRED)
  add_definitions( -DMPICH_IGNORE_CXX_SEEK )
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

#INCLUDE(CTest)

#find_library(DL_LIB dl)

set (SURFLIB_SRCS src/Rocsurf.C src/Manifold_2.C src/Generic_element_2.C src/interpolate_to_centers.C 
                 src/compute_element_normals.C src/compute_element_areas.C src/compute_bounded_volumes.C src/compute_curvature.C)
set (TEST_SRCS test/surfmap.C test/surftest1.C test/surftest2.C test/imtest.C test/im2hdf.C test/smf2hdf.C test/neu2hdf.C
     test/mcntest.C test/obj2hdf.C test/msh2hdf.C test/serialize.C)
set (ALL_SURF_SRCS "${SURFLIB_SRCS} ${TEST_SRCS}")

#set (TEST_SRCS test/surftest.C)
set_source_files_properties(${ALL_SURF_SRCS} PROPERTIES COMPILE_FLAGS "-fPIC ${MPI_CXX_COMPILE_FLAGS}" )
#set(TEST_SRCS src/COMTest.C)

# rpath settings
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include_directories(include test)

IF(NOT BUILD_STATIC)
  add_library(Rocsurf SHARED ${SURFLIB_SRCS})
ELSE()
  add_library(Rocsurf STATIC ${LIB_SRCS})
  add_definitions( -DSTATIC_LINK )
ENDIF()
target_link_libraries(Rocsurf Rocmap Rocblas Rocout ${MPI_CXX_LIBRARIES})

add_executable(surfmap test/surfmap.C)
target_link_libraries(surfmap Rocsurf ${MPI_CXX_LIBRARIES} ${HDF5_LIB} )
add_executable(surftest1 test/surftest1.C)
target_link_libraries(surftest1 Rocsurf ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
add_executable(surftest2 test/surftest2.C)
target_link_libraries(surftest2 Rocsurf ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
add_executable(imtest test/imtest.C)
target_link_libraries(imtest Rocsurf ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
add_executable(im2hdf test/im2hdf.C)
target_link_libraries(im2hdf Rocsurf ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
add_executable(smf2hdf test/smf2hdf.C)
target_link_libraries(smf2hdf Rocsurf ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
add_executable(neu2hdf test/neu2hdf.C)
target_link_libraries(neu2hdf Rocsurf ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
add_executable(mcntest test/mcntest.C)
target_link_libraries(mcntest Rocsurf ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
add_executable(obj2hdf test/obj2hdf.C)
target_link_libraries(obj2hdf Rocsurf ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
add_executable(msh2hdf test/msh2hdf.C)
target_link_libraries(msh2hdf Rocsurf ${MPI_CXX_LIBRARIES} ${HDF5_LIB})
add_executable(serialize test/serialize.C)
target_link_libraries(serialize Rocsurf ${MPI_CXX_LIBRARIES} ${HDF5_LIB})

IF(MPI_LINK_FLAGS)
  SET_TARGET_PROPERTIES(Rocsurf surfmap surftest1 surftest2 imtest im2hdf
    smf2hdf neu2hdf mcntest obj2hdf msh2hdf serialize PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
ENDIF()

#ADD_TEST(RunAllTests ${EXECUTABLE_OUTPUT_PATH}/testx testresults.txt)
#ADD_TEST(TestObject:ValidHandle ${TEST_RESULTS} TestObject:ValidHandle testresults.txt)
#ADD_TEST(TestObject:Exists ${TEST_RESULTS} TestObject:Exists testresults.txt)

INSTALL(TARGETS Rocsurf RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
