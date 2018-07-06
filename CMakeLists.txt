cmake_minimum_required (VERSION 2.8)
project(Rocstar CXX C Fortran)
set(BUILD_STATIC FALSE CACHE BOOL "Build static libraries")
set(ENABLE_MPI TRUE CACHE BOOL "Build with MPI Support")
set(ENABLE_CGNS TRUE CACHE BOOL "Build with CGNS Support")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single directory for all executables.")
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Single directory for all libraries and archives.")
mark_as_advanced (LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

if(NOT ENABLE_MPI)
  add_definitions( -DDUMMY_MPI )
else()
  find_package(MPI REQUIRED)
  add_definitions( -DMPICH_IGNORE_CXX_SEEK )
  include_directories(${MPI_INCLUDE_PATH})
  set(FC mpif90)
endif()

if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

include(CTest)
if(BUILD_TESTING)
  set(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  mark_as_advanced(BUILDNAME)
endif(BUILD_TESTING)

# finding IMPACT libs and incs
find_library(SITCOM_LIB SITCOM)
get_filename_component(IMPACT_LIBPATH ${SITCOM_LIB} PATH)
set(IMPACT_BINPATH "${IMPACT_LIBPATH}/../bin" )
find_library(SITCOMF_LIB SITCOMF HINTS ${IMPACT_LIBPATH})
find_library(SIMIN_LIB SimIN HINTS ${IMPACT_LIBPATH})
find_library(SIMOUT_LIB SimOUT HINTS ${IMPACT_LIBPATH})
find_library(SIMPAL_LIB Simpal HINTS ${IMPACT_LIBPATH})
find_library(SIM_LIB SIM HINTS ${IMPACT_LIBPATH})
find_library(SURFMAP_LIB SurfMap HINTS ${IMPACT_LIBPATH})
find_library(SURFUTIL_LIB SurfUtil HINTS ${IMPACT_LIBPATH})
find_library(SURFX_LIB SurfX HINTS ${IMPACT_LIBPATH})
set(IMPACT_LIB ${SITCOM_LIB} ${SITCOMF_LIB} ${SIMIN_LIB}
               ${SIMIN_LIB} ${SIMOUT_LIB} ${SIMPAL_LIB}
               ${SIM_LIB} ${SURFMAP_LIB} ${SURFUTIL_LIB}
               ${SURFX_LIB})
mark_as_advanced(IMPACT_LIB SITCOMF_LIB SIMIN_LIB SIMOUT_LIB SIMPAL_LIB SIM_LIB SURFMAP_LIB SURFUTIL_LIB SURFX_LIB)
find_file(IMPACT_HDR com.h)
get_filename_component(IMPACT_INCPATH ${IMPACT_HDR} PATH)
set(IMPACT_INC ${IMPACT_INCPATH})
set(ROCCOM_INC ${IMPACT_INCPATH})
set(ROCIN_INC ${IMPACT_INCPATH})
set(SIMPAL_INC ${IMPACT_INCPATH})
set(SURF_INC ${IMPACT_INCPATH})
set(MAP_INC ${IMPACT_INCPATH})
include_directories(${IMPACT_INCPATH})
mark_as_advanced(IMPACT_INC ROCCOM_INC ROCIN_INC SIMPAL_INC SURF_INC MAP_INC)

# other needed libraries
find_library(HDF4_LIB NAMES hdf4 libdf df PATHS /usr/lib /usr/lib64 HINTS df libdf hdf)
find_library(MFHDF_LIB NAMES mfhdf4 mfhdf PATHS /usr/lib /usr/lib64 HINTS mfhdf mfhdf4)
find_path(HDF4_INC hdf.h HINTS /usr/include /usr/include/hdf)

if (ENABLE_CGNS)
  find_package(HDF5 REQUIRED) # will set HDF5_LIBRARIES
  find_library(CGNS_LIB cgns HINTS /usr/lib/ /usr/lib64/)
  add_definitions(-DUSE_CGNS)
endif()

# compiler definitions
if(MPI_CXX_COMPILER)
  add_definitions(-D_ROCSTAR_PARALLEL_)
endif()

# original rocstar stuff below here
add_subdirectory(Rocprop)
add_subdirectory(Rocmop)
add_subdirectory(Rocon)
add_subdirectory(Rocman)
add_subdirectory(RocfluidMP)
add_subdirectory(Rocfrac)
add_subdirectory(Rocburn)
add_subdirectory(Rocprep)
add_subdirectory(testing)

include_directories(./Roccom/include ./Rocman/include)
add_executable(rocstar src/rocstar.C)
set_source_files_properties(src/rocstar.C PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
set_target_properties(rocstar PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
target_link_libraries(rocstar Rocman ${IMPACT_LIB} ${MPI_CXX_LIBRARIES})

install(TARGETS rocstar RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

# Add target to build the documentation
option(BUILD_DOCUMENTATION "Create and install the HTML based API
documentation (requires Doxygen)" ON)
if(BUILD_DOCUMENTATION)
  find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation.")
  endif()
  set( doxyfile_in          ${CMAKE_CURRENT_SOURCE_DIR}/Docs/Rocstar.cfg.in )
  set( doxyfile             ${CMAKE_CURRENT_BINARY_DIR}/Rocstar.cfg        )
  set( documentation        ${CMAKE_CURRENT_BINARY_DIR}/Docs                 )
  configure_file( ${doxyfile_in} ${doxyfile} @ONLY )
  add_custom_command( OUTPUT ${documentation}
                      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
                      MAIN_DEPENDENCY ${doxyfile} ${doxyfile_in}
                      COMMENT "Generating documentation")
  add_custom_target( Docs DEPENDS ${documentation} )
  install( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Docs DESTINATION share/Rocstar )
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Docs/RocstarQuickStart.pdf DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Docs/)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Docs/legacy DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Docs/)
endif()

include(CPack)
