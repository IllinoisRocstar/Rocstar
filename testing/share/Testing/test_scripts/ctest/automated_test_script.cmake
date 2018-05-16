# ====== PLATFORM SPECIFIC SETTINGS =======
# THESE VARIABLES MUST BE SET FOR EACH PLATFORM ON WHICH TESTING IS ENABLED
# The following variables are inherited from the environment
# --------
# SVNCOMMAND = "/path/to/svn" 
# PROJECT_SOURCE = "/path/to/project/source" (e.g. /home/mtcampbe/Nightly/Rocstar)
# PROJECT_ROOT = "repository/path/to/project" (e.g. Rocstar/trunk)
# These are picked up from the user's environment
# CMAKECOMMAND = "/home/mtcampbe/VTK/Install/bin/cmake"
# CTESTCOMMAND = "/home/mtcampbe/VTK/Install/bin/ctest"
# This setting is optional
# PROJECT_CONFIGURATION_OPTIONS = "-DBOOSTROOT=/path/to/boost"
# --------
SET (_IR_REPO_ "file:///Projects/IR/SourceRepository")
SET (_PATH_TO_SVN_ $ENV{SVNCOMMAND})
SET (_PATH_TO_PROJECT_SOURCE_ $ENV{PROJECT_SOURCE})
SET (_PROJECT_ROOT_ $ENV{PROJECT_ROOT})
SET (_CTEST_COMMAND_ $ENV{CTESTCOMMAND})
SET (_CMAKE_COMMAND_ $ENV{CMAKECOMMAND})
SET (_PATH_TO_PROJECT_BINARY_ $ENV{PROJECT_BUILD_DIRECTORY})
# OPTIONAL
SET (_CMAKE_CONFIGURATION_OPTIONS_ ${PROJECT_CONFIGURATION_OPTIONS})
# =========================================

# -----------------------------------------------------------  
# -- Get environment
# -----------------------------------------------------------  

## -- Set hostname
## --------------------------
find_program(HOSTNAME_CMD NAMES hostname)
macro(gethostname name flag)
  exec_program("${HOSTNAME_CMD}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(gethostname)

gethostname(shorthost -s)
exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)
#exec_program(${HOSTNAME_CMD} -s OUTPUT_VARIABLE shorthost)
set(CTEST_SITE                          "${HOSTNAME}")
## -- Set site / build name
## --------------------------

find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)

getuname(osname -s)
getuname(osrel  -r)
getuname(cpu    -m)

set(CTEST_BUILD_NAME                    "${osname}-${cpu}")

SET (CTEST_SOURCE_DIRECTORY "${_PATH_TO_PROJECT_SOURCE_}")
SET (CTEST_BINARY_DIRECTORY "${_PATH_TO_PROJECT_BINARY_}")
# should ctest wipe the binary tree before running
SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

SET (CTEST_SVN_COMMAND "${_PATH_TO_SVN_}")
SET (CTEST_SVN_CHECKOUT  "${CTEST_SVN_COMMAND} co ${_IR_REPO_}/${_PROJECT_ROOT_} ${CTEST_SOURCE_DIRECTORY}")
# set any extra environment variables to use during the execution of the script here:
#SET (CTEST_ENVIRONMENT
# "CMAKE_PREFIX_PATH=${_CMAKE_PREFIX_PATH_}"
# "BOOST_ROOT=${_PATH_TO_BOOST_}"
#)

# which ctest command to use for running the dashboard
SET(MODEL Experimental)
IF(${CTEST_SCRIPT_ARG} MATCHES Nightly)
  SET(MODEL Nightly)
ENDIF(${CTEST_SCRIPT_ARG} MATCHES Nightly)
IF(${CTEST_SCRIPT_ARG} MATCHES Continuous)
  SET(MODEL Continuous)
ENDIF(${CTEST_SCRIPT_ARG} MATCHES Continuous)
SET (CTEST_COMMAND "${_CTEST_COMMAND_} -D ${MODEL}")


find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
#SET(IR_SOURCE_REPOSITORY "file:///Projects/IR/SourceRepository")
if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  set(CTEST_CHECKOUT_COMMAND "${CTEST_SVN_COMMAND} co ${_IR_REPO_}/${_PROJECT_ROOT_} ${CTEST_SOURCE_DIRECTORY}")
endif()

# what cmake command to use for configuring this dashboard
SET (CTEST_CMAKE_COMMAND "${_CMAKE_COMMAND_}")


####################################################################
# The values in this section are optional you can either
# have them or leave them commented out
####################################################################

# this is the initial cache to use for the binary tree, be careful to escape
# any quotes inside of this string if you use it
SET (CTEST_INITIAL_CACHE "
CMAKE_GENERATOR:INTERNAL=Unix Makefiles
BUILDNAME:STRING=${CTEST_BUILD_NAME}
SITE:STRING=${CTEST_SITE}
SVN_UPDATE_OPTIONS:STRING=update
")

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_UPDATE_COMMAND "${CTEST_SVN_COMMAND}")
set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} ${_CMAKE_CONFIGURATION_OPTIONS_} ${CTEST_SOURCE_DIRECTORY}")
#set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} -DCMAKE_PREFIX_PATH=/Projects/Tomography_code/Support/Mercury -DBoost_NO_SYSTEM_PATHS=ON -DBoost_NO_BOOST_CMAKE=ON ${CTEST_SOURCE_DIRECTORY}")
#set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DWITH_TESTING:BOOL=ON ${CTEST_BUILD_OPTIONS}")
#set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"-G${CTEST_CMAKE_GENERATOR}\"")
#set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"${CTEST_SOURCE_DIRECTORY}\"")

ctest_start("${MODEL}")
ctest_update()
ctest_configure()
ctest_build()
ctest_test()
if (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
  ctest_coverage()
endif (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
if (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
  ctest_memcheck()
endif (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
ctest_submit()
