cmake_minimum_required (VERSION 2.8)
project(Rocface)

load_cache(../)
IF(NOT ENABLE_MPI)
  add_definitions( -DDUMMY_MPI )
ELSE()
  FIND_PACKAGE(MPI REQUIRED)
  add_definitions( -DMPICH_IGNORE_CXX_SEEK )
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

#INCLUDE(CTest)
find_path(ROCCOM_INC roccom.h HINTS ../Roccom/include) 
find_path(ROCIN_INC HDF4.h HINTS ../Roccom/Rocin/include) 

set (RFC_SRCS src/Base/RFC_Window_base.C src/Base/RFC_Window_base_IO.C src/Base/RFC_Window_base_IO_tecplot.C 
              src/Base/RFC_Window_base_IO_binary.C src/Overlay/RFC_Window_overlay.C src/Overlay/RFC_Window_overlay_fea.C
              src/Overlay/RFC_Window_overlay_IO.C src/Overlay/Overlay_primitives.C src/Overlay/Overlay_0d.C 
              src/Overlay/Overlay_init.C src/Overlay/Overlay.C src/Overlay/Overlay_IO.C src/Base/Vector_n.C
              src/Transfer/RFC_Window_transfer.C src/Transfer/RFC_Window_transfer_comm.C src/Transfer/Transfer_base.C
              src/Transfer/Transfer_2f.C src/Transfer/Transfer_2n.C src/Rocface.C)
#              src/Base/writer.C src/Base/rfc_assertions.C)
#              src/Transfer/Transfer_2f.C src/Transfer/Transfer_2n.C src/Rocface.C src/Overlay/Triangulation.C 
set (TEST_SRCS test/ex1.C test/ex2.C test/ex3.C test/ex4.C test/ex5.C test/plot3d.C test/reptrans.C test/readsdv.C test/ex6DataTransfer.C test/meshio.C)
set (UTIL_SRCS util/surfdiver.C util/rfctest.C util/autosurfer.C util/surfjumper.C util/surfextractor.C)
set (CGAL_SRCS CGAL/src/Color.C CGAL/src/Double.C CGAL/src/Origin.C CGAL/src/aff_transformation_tags.C CGAL/src/assertions.C CGAL/src/io.C)
set (LIB_SRCS ${RFC_SRCS} ${CGAL_SRCS})
set (ALL_SRCS ${LIB_SRCS} ${TEST_SRCS} ${UTIL_SRCS})

set_source_files_properties(${ALL_SRCS} PROPERTIES COMPILE_FLAGS "-fPIC ${MPI_CXX_COMPILE_FLAGS}" )
#set(TEST_SRCS src/MANTest.C)

# rpath settings
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include_directories(include ${ROCCOM_INC} ${ROCIN_INC} include/Base include/Transfer include/Overlay CGAL/include)


IF(NOT BUILD_STATIC)
  add_library(Rocface SHARED ${LIB_SRCS})
ELSE()
  add_library(Rocface STATIC ${LIB_SRCS})
  add_definitions( -DSTATIC_LINK )
ENDIF()
target_link_libraries(Rocface Rocsurf Rocmap Rocsurf Roccom ${MPI_CXX_LIBRARIES})

add_executable(ex1 test/ex1.C test/meshio.C)
add_executable(ex2 test/ex2.C test/meshio.C)
add_executable(ex3 test/ex3.C test/meshio.C)
add_executable(ex4 test/ex4.C test/meshio.C)
add_executable(ex5 test/ex5.C test/meshio.C)
add_executable(ex6DataTransfer test/ex6DataTransfer.C test/meshio.C)
add_executable(plot3d test/plot3d.C test/meshio.C)
add_executable(reptrans test/reptrans.C test/meshio.C)
add_executable(readsdv test/readsdv.C test/meshio.C)
target_link_libraries(ex1 Rocface ${HDF5_LIB})
target_link_libraries(ex2 Rocface ${HDF5_LIB})
target_link_libraries(ex3 Rocface ${HDF5_LIB})
target_link_libraries(ex4 Rocface ${HDF5_LIB})
target_link_libraries(ex5 Rocface ${HDF5_LIB})
target_link_libraries(ex6DataTransfer Rocface ${HDF5_LIB})
target_link_libraries(plot3d Rocface ${HDF5_LIB})
target_link_libraries(reptrans Rocface ${HDF5_LIB})
target_link_libraries(readsdv Rocface ${HDF5_LIB})

add_executable(surfdiver util/surfdiver.C)
target_link_libraries(surfdiver Rocface ${HDF5_LIB})
add_executable(rfctest util/rfctest.C)
target_link_libraries(rfctest Rocface ${HDF5_LIB})
add_executable(autosurfer util/autosurfer.C)
target_link_libraries(autosurfer Rocface ${HDF5_LIB})
add_executable(surfjumper util/surfjumper.C)
target_link_libraries(surfjumper Rocface ${HDF5_LIB})
add_executable(surfextractor util/surfextractor.C)
target_link_libraries(surfextractor Rocface ${HDF5_LIB})

SET_TARGET_PROPERTIES(Rocface PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(ex1 PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(ex2 PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(ex3 PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(ex4 PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(ex5 PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(ex6DataTransfer PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(plot3d PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(reptrans PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(readsdv PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(surfdiver PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(rfctest PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(autosurfer PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(surfjumper PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
SET_TARGET_PROPERTIES(surfextractor PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

#ADD_TEST(RunAllTests ${EXECUTABLE_OUTPUT_PATH}/testx testresults.txt)
#ADD_TEST(TestObject:ValidHandle ${TEST_RESULTS} TestObject:ValidHandle testresults.txt)
#ADD_TEST(TestObject:Exists ${TEST_RESULTS} TestObject:Exists testresults.txt)

INSTALL(TARGETS Rocface surfdiver rfctest autosurfer surfjumper surfextractor RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
