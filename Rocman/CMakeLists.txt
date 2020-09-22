cmake_minimum_required(VERSION 3.0)

if(NOT BUILD_SHARED_LIBS)
  add_definitions(-DSTATIC_LINK)
endif()

# rpath settings
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
add_library(Rocman
    src/rocstar_driver.C
    # Couplings
    src/RocstarCoupling.C
    src/builtin_couplings.C
    src/derived_couplings.C
    # Actions
    src/basic_actions.C
    src/SurfDiver.C
    src/transfer_actions.C
    # Agents
    src/RocstarAgent.C
    src/SolidAgent.C
    src/BurnAgent.C
    src/FluidAgent.C
    # Schedulers
    # Utility
    src/Control_parameters.C
    src/rocman.C
)

target_link_libraries(Rocman IMPACT::SIM IMPACT::SITCOM IMPACT::SITCOMF)

target_include_directories(Rocman
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/rocstar>
)
set_target_properties(Rocman PROPERTIES VERSION ${ROCSTAR_VERSION}
        SOVERSION ${ROCSTAR_MAJOR_VERSION})

install(DIRECTORY include/ 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/rocstar)
install(TARGETS Rocman
    EXPORT Rocstar
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
