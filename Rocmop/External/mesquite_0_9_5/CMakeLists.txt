cmake_minimum_required (VERSION 2.8)
project(Mesquite)

FILE(GLOB LIB_SRCS srcLinks/*.cpp)

# rpath settings
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include_directories(include includeLinks)

add_definitions(-DMESQUITE)
add_definitions(-DUSE_STD_INCLUDES -DUSE_C_PREFIX_INCLUDES -DHAVE_VSNPRINTF -DHAVE_VSPRINTF -DHAVE_CLOCK -DHAVE_TIMES -DHAVE_CBRT)

IF(NOT BUILD_STATIC)
  add_library(mesquite_95 SHARED ${LIB_SRCS})
ELSE()
  add_library(mesquite_95 STATIC ${LIB_SRCS})
  add_definitions( -DSTATIC_LINK )
ENDIF()

INSTALL(TARGETS mesquite_95 RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
