cmake_minimum_required(VERSION 3.0)

set(LIB_SRCS
    src/Propagation_3.C
    src/MarkerParticles_3.C
    src/FaceOffset_3.C
    src/Rocprop.C
    src/NuLaplacian.C
    src/cons_diff.C
    src/AnisotropicSmoothing.C
    src/quadric_analysis.C
    src/detect_features.C)
#set(TEST_SRCS
#    test/proptest.C
#    test/mptest.C
#    test/compare_normals.C
#    test/rotate.C
#    test/advectest.C
#    test/detectfeas.C
#    test/PointPropagate.C)
#set(TEST_SRCS src/MANTest.C)

add_library(Rocprop ${LIB_SRCS})
target_include_directories(Rocprop 
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/rocstar>)
target_link_libraries(Rocprop IMPACT::SurfUtil)
set_target_properties(Rocprop PROPERTIES VERSION ${ROCSTAR_VERSION}
        SOVERSION ${ROCSTAR_MAJOR_VERSION})

#add_executable(proptest test/proptest.C)
#target_link_libraries(proptest Rocprop)
#add_executable(mptest test/mptest.C)
#target_link_libraries(mptest Rocprop)
#add_executable(compare_normals test/compare_normals.C)
#target_link_libraries(compare_normals Rocprop)
#add_executable(rotate test/rotate.C)
#target_link_libraries(rotate Rocprop)
#add_executable(advectest test/advectest.C test/PointPropagate.C)
#target_link_libraries(advectest Rocprop)
#add_executable(detectfeas test/detectfeas.C)
#target_link_libraries(detectfeas Rocprop)

install(TARGETS Rocprop
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
