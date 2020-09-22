cmake_minimum_required(VERSION 3.0)

# set project metadata #########################################################

set(ROCSTAR_MAJOR_VERSION 5)
set(ROCSTAR_MINOR_VERSION 0)
set(ROCSTAR_PATCH_VERSION 1)
set(ROCSTAR_EXTRA_VERSION "" CACHE STRING "ROCSTAR extra version string")
set(ROCSTAR_VERSION "${ROCSTAR_MAJOR_VERSION}.${ROCSTAR_MINOR_VERSION}")
set(ROCSTAR_VERSION
    "${ROCSTAR_VERSION}.${ROCSTAR_PATCH_VERSION}${ROCSTAR_EXTRA_VERSION}")
set(ROCSTAR_SHORT_LICENSE "University of Illinois/NCSA Open Source License")
set(ROCSTAR_HOMEPAGE_URL "http://git.illinois.rocstar/rocstar_modern/rocstar")
set(ROCSTAR_AUTHOR
    "Masoud Safdari"
    "Alessandro Gondolo"
    "Michael Pope"
    "Eric Shaffer"
    "Fady Najjar"
    "Andreas Haselbacher"
    "Xiangmin Jiao"
    "Michael Campbell")
set(ROCSTAR_AUTHOR_EMAIL "msafdari@illinoisrocstar.com")
set(ROCSTAR_DESCRIPTION "Rocstar is a multiphysics simulation application \
designed for coupled multiphysics simulations involving fluid-structure \
interaction (FSI) across moving, reacting interfaces. Rocstar couples multiple \
domain-specific simulation packages and disparately discretized domains and \
provides several simulation-supporting services including conservative and \
accurate data transfer, surface propagation, and parallel I/O. Rocstar is \
MPI-parallel. Rocstar was originally developed at the University of Illinois \
Center for Simulation of Advanced Rockets (CSAR) under Department of Energy \
ASCI program funding.")
set(PROJECT_NAME "Rocstar version ${ROCSTAR_VERSION}")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set project defaults #########################################################

project(${PROJECT_NAME} CXX C Fortran)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single
 directory for all executables.")
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Single directory
 for all libraries and archives.")
mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build.")
endif()
# Set the possible values of build type for cmake-gui
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")

# Installation Directories
include(GNUInstallDirs)

# rpath settings
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Set configure options ########################################################

option(BUILD_SHARED_LIBS "Build shared libraries." ON)
option(ENABLE_TESTING "Build with tests." OFF)
option(ENABLE_PROFILING "Build with profiling Support" OFF)
option(ENABLE_DEBUG "Enable debug print statements" OFF)

#option(ENABLE_MPI "Build with MPI support." ON) # Set to always ON
set(ENABLE_MPI ON)

# Set warnings #################################################################

# Turn on default warnings for GNU compilers.
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Basic warnings for compatibility and conformance to standard
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic")
  set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -std=c11 -Wpedantic")
  set(CMAKE_Fortran_FLAGS "${CMAKE_CXX_FLAGS} -std=f2008 -Wpedantic")

  # Additional warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual -Wnon-virtual-dtor")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wstrict-null-sentinel -Wredundant-decls")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfloat-conversion")
  endif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9)

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wno-unused")

#[[
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -Wno-unused")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 6) # Bugged till gcc-6
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-dummy-argument")
  endif()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wcharacter-truncation")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wshadow")
]]

  # Fortran optimization warnings.
  #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Warray-temporaries")
endif()

# Find libraries ##########$####################################################

# Finding IMPACT libs
find_package(IMPACT REQUIRED 2.1.3)

# Find MPI package or use a dummy.
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_INCLUDE_PATH})
  if(MPI_LIBRARY)
    set(MPI_LIBRARIES ${MPI_LIBRARY})
  endif(MPI_LIBRARY)
  if(MPI_EXTRA_LIBRARY)
    set(MPI_LIBRARIES ${MPI_LIBRARIES} "${MPI_EXTRA_LIBRARY}")
  endif(MPI_EXTRA_LIBRARY)
  # Needed for mpich 2
  add_definitions(-DMPICH_IGNORE_CXX_SEEK)
else()
  add_definitions(-DDUMMY_MPI)
endif()

if(ENABLE_DEBUG)
  add_definitions(-DENABLE_DEBUG)
  # check_datastruct flag does more checking inside rocstar and prints out
  # debug information
  add_definitions(-DCHECK_DATASTRUCT)
  # check gradients within RFLU/RFLO
  add_definitions(-DCHECK_GRAD)
  # rocpart debug & CE debug information
  if(ENABLE_MPI)
    add_definitions(-DPLAG_MPI_DEBUG -DPLAG_CECELLS_MPI_DEBUG)
  else(ENABLE_MPI)
    add_definitions(-DPLAG_CECELLS_DEBUG -DPLAG_DEBUG)
  endif(ENABLE_MPI)
endif(ENABLE_DEBUG)

if(ENABLE_PROFILING)
  find_library(GPROF_LIB tcmalloc_and_profiler)
  get_filename_component(GPROF_LNKPATH ${GPROF_LIB} PATH)
  link_directories(${GPROF_LNKPATH})
  find_file(GPROF_HDR heap-profiler.h path /usr/include/gperftools/)
  get_filename_component(GPROF_INCPATH ${GPROF_HDR} PATH)
  include_directories(AFTER ${GPROF_INCPATH})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_GPROF")
  # ROCPROF is an internal profiler which can be turned on if desired
  # add_definitions(-DROCPROF)
endif(ENABLE_PROFILING)

# Add project submodules #######################################################

# compiler definitions/flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

# Original rocstar stuff below here
add_subdirectory(Rocman)
add_subdirectory(Rocmop)
add_subdirectory(Rocburn)
add_subdirectory(RocfluidMP)
add_subdirectory(Rocprep)
add_subdirectory(Rocprop)

# Perform Testing ##############################################################

if(ENABLE_TESTING)
  include(CTest)
  enable_testing()
  add_subdirectory(testing)
  set(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard.")
  mark_as_advanced(BUILDNAME)
endif(ENABLE_TESTING)

# Build rocstar executable #####################################################

add_executable(rocstar src/rocstar.C)
target_link_libraries(rocstar Rocman ${GPROF_LIB})

# Installation #################################################################

install(TARGETS rocstar
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Documentation ################################################################

# Add target to build the documentation
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation
 (requires Doxygen)" ON)
if(BUILD_DOCUMENTATION)
  find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation.")
  endif()
  set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Docs/Rocstar.cfg.in)
  set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Rocstar.cfg)
  set(documentation ${CMAKE_CURRENT_BINARY_DIR}/Docs)
  configure_file(${doxyfile_in} ${doxyfile} @ONLY)
  add_custom_command(OUTPUT ${documentation}
      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
      MAIN_DEPENDENCY ${doxyfile} ${doxyfile_in}
      COMMENT "Generating documentation")
  add_custom_target(Docs DEPENDS ${documentation})
endif(BUILD_DOCUMENTATION)

# Install the Docs directory regardless of building the doxygen files
# Rocstar comes with documentation from the Git repository already
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Docs/legacy
    ${CMAKE_CURRENT_SOURCE_DIR}/Docs/RocstarQuickStart.pdf
    ${CMAKE_CURRENT_SOURCE_DIR}/Docs/RocstarTutorial.pdf
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Docs/)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Docs DESTINATION ${CMAKE_INSTALL_DOCDIR})

# Packaging ####################################################################

set(CPACK_PACKAGE_NAME "Rocstar")
set(CPACK_PACKAGE_VENDOR "${ROCSTAR_AUTHOR}")
set(CPACK_PACKAGE_VERSION "${ROCSTAR_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${ROCSTAR_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR "${ROCSTAR_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH "${ROCSTAR_PATCH_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION "${ROCSTAR_DESCRIPTION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${ROCSTAR_DESCRIPTION}")
set(CPACK_PACKAGE_HOMEPAGE_URL "")
set(CPACK_PACKAGE_CONTACT "${ROCSTAR_AUTHOR_EMAIL}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")

file(READ "/etc/os-release" OS_INFO)
string(REGEX MATCH "Ubuntu|CentOS" DIST ${OS_INFO})
if (DIST STREQUAL "Ubuntu")
  set(Debian TRUE)
elseif (DIST STREQUAL "CentOS")
  set(CentOS TRUE)
endif ()

if(WIN32)
  set(CPACK_GENERATOR NSIS)
elseif (Debian)
  if(CMAKE_VERSION VERSION_LESS 3.10)
    message(WARNING "CMake 3.10 or above is required to build a DEB package.")
  else()
    set(CPACK_GENERATOR DEB)

    # set Debian Packaging settings
    #set(CPACK_DEBIAN_PACKAGE_NAME "")
    #set(CPACK_DEBIAN_PACKAGE_VERSION "")
    #set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "")
    #set(CPACK_DEBIAN_PACKAGE_MAINTAINER "")
    #set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "")
    set(CPACK_DEBIAN_PACKAGE_SECTION "science")
    #set(CPACK_DEBIAN_COMPRESSION_TYPE "gzip")
    #set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
    set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "${CPACK_PACKAGE_HOMEPAGE_URL}")
    set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
    set(CPACK_DEBIAN_PACKAGE_DEBUG ON)  # For testing only
    #set(CPACK_DEBIAN_PACKAGE_PREDEPENDS "")
    #set(CPACK_DEBIAN_PACKAGE_ENHANCES "")
    #set(CPACK_DEBIAN_PACKAGE_BREAKS "")
    #set(CPACK_DEBIAN_PACKAGE_CONFLICTS "")
    #set(CPACK_DEBIAN_PACKAGE_PROVIDES "")
    #set(CPACK_DEBIAN_PACKAGE_REPLACES "")
    #set(CPACK_DEBIAN_PACKAGE_RECOMMENDS "")
    #set(CPACK_DEBIAN_PACKAGE_SUGGESTS "")
    #set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "")

    set(CPACK_DEBIAN_PACKAGE_CONTROL_STRICT_PERMISSION TRUE)  # CMake 3.4
    #set(CPACK_DEBIAN_PACKAGE_SOURCE "")  # CMake 3.5
    #set(CPACK_DEB_COMPONENT_INSTALL OFF)  # CMake 3.6
    set(CPACK_DEBIAN_FILE_NAME "DEB-DEFAULT")  # CMake 3.6
    #set(CPACK_DEBIAN_PACKAGE_RELEASE "1")  # CMake 3.6
    set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)  # CMake 3.6
    #set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS_POLICY "=")  # CMake 3.6
    #set(CPACK_DEBIAN_ARCHIVE_TYPE "gnutar")  # CMake 3.7
    #set(CPACK_DEBIAN_PACKAGE_EPOCH "")  # Cmake 3.10
  endif()
elseif (CentOS)
  set(CPACK_GENERATOR RPM)

  # set RPM package settings
  # set(CPACK_RPM_COMPONENT_INSTALL "") for multiple packages
  # set(CPACK_RPM_PACKAGE_SUMMARY "")
  # set(CPACK_RPM_PACKAGE_NAME "")
  set(CPACK_RPM_FILE_NAME "RPM-DEFAULT")
  # set(CPACK_RPM_MAIN_COMPONENT "")
  # set(CPACK_RPM_PACKAGE_EPOCH "") for package numbering when the version # doesn't change
  # set(CPACK_RPM_PACKAGE_VERSION "")
  # set(CPACK_RPM_PACKAGE_ARCHITECTURE "")
  # set(CPACK_RPM_PACKAGE_RELEASE "")
  # set(CPACK_RPM_PACKAGE_RELEASE_DIST "")
  # set(CPACK_RPM_PACKAGE_LICENSE "") # By default it sets this to "unknown"
  set(CPACK_RPM_PACKAGE_GROUP "science")
  set(CPACK_RPM_PACKAGE_VENDOR "IllinoisRocstar")
  set(CPACK_RPM_PACKAGE_URL "${CPACK_PACKAGE_HOMEPAGE_URL}")
  # set(CPACK_RPM_PACKAGE_DESCRIPTION "")
  # set(CPACK_RPM_COMPRESSION_TYPE "")
  set(CPACK_RPM_PACKAGE_AUTOREQ ON) # Enable automatic shared lib dep detection
  set(CPACK_RPM_PACKAGE_AUTOPROV ON) # Enable automatic listing of shared libs
  set(CPACK_RPM_PACKAGE_AUTOREQPROV ON)
  # set(CPACK_RPM_PACKAGE_REQUIRES "")
  # set(CPACK_RPM_PACKAGE_CONFLICTS "")
  # set(CPACK_RPM_PACKAGE_REQUIRES_PRE "")
  # set(CPACK_RPM_PACKAGE_REQUIRES_POST "")
  # set(CPACK_RPM_PACKAGE_REQUIRES_POSTUN "")
  # set(CPACK_RPM_PACKAGE_REQUIES_PREUN "")
  # set(CPACK_RPM_PACKAGE_SUGGESTS "")
  # set(CPACK_RPM_PACKAGE_PROVIDES "")
  # set(CPACK_RPM_PACKAGE_OBSOLETES "")
  # set(CPACK_RPM_PACKAGE_RELOCATABLE OFF)
  # set(CPACK_RPM_SPEC_INSTALL_POST "")
  # set(CPACK_RPM_SPEC_MORE_DEFINE "")
  set(CPACK_RPM_PACKAGE_DEBUG ON)
  # set(CPACK_RPM_USER_BINARY_SPECFILE "")
  # set(CPACK_RPM_GENERATE_USER_BINARY_SPECFILE_TEMPLATE "")
  # set(CPACK_RPM_PRE_INSTALL_SCRIPT_FILE "")
  # set(CPACK_RPM_POST_INSTALL_SCRIPT_FILE "")
  # set(CPACK_RPM_USER_FILELIST "")
  # set(CPACK_RPM_CHANGELOG_FILE "")
  # set(CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST "")
  # set(CPACK_RPM_DEFAULT_USER "") User & group should always be root
  # set(CPACK_RPM_DEFAULT_GROUP "")
else ()
    message(WARNING "Package generation is currently only supported on CentOS and Ubuntu distributions")
endif ()
# Lists exe and text label used for Start Menu shortcuts
#set(CPACK_PACKAGE_EXECUTABLES "")
# Must have matching Start Menu shortcut to create a Desktop link
#set(CPACK_CREATE_DESKTOP_LINKS "")

# WIN32 NOTE: - Run CPack from a MSVC environment as dumpbin is required
#               for the fixup_bundle search of DLL dependecies.
#             - The system libraries pulled in are release DLLs.
#               Do NOT package in Debug! THAT APPLIES TO DEPENDENCIES TOO!!
include(CPack)
