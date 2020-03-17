cmake_minimum_required(VERSION 3.0)
project(Mesquite)

file(GLOB LIB_SRCS srcLinks/*.cpp)

add_definitions(
#    -DUSE_STD_INCLUDES
#    -DUSE_C_PREFIX_INCLUDES
    -DHAVE_VSNPRINTF
    -DHAVE_VSPRINTF
    -DHAVE_CLOCK
    -DHAVE_TIMES
    -DHAVE_CBRT)

add_library(mesquite_95 ${LIB_SRCS})
target_include_directories(mesquite_95 SYSTEM PUBLIC include includeLinks)

install(TARGETS mesquite_95
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
