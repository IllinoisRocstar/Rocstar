cmake_minimum_required(VERSION 3.0)
# project(Rocburn C CXX Fortran)

set(CORE_SRCS
    01.burn_interface_data.f90
    calcdist.f90
    rocburn_2D.f90
    Rocburn_load_module.f90)
set(APN_SRCS
    Rocburn_APN/02.burn_APN_global_data.f90
    Rocburn_APN/APN_input_0d.f90
    Rocburn_APN/Rocburn_1D_APN.f90)
set(ZN_SRCS
    Rocburn_ZN/02.burn_ZN_global_data.f90
    Rocburn_ZN/Rocburn_1D_ZN.f90
    Rocburn_ZN/ZN_calc_burning_rate.f90
    Rocburn_ZN/ZN_gen_grid.f90
    Rocburn_ZN/ZN_input_0d.f90
    Rocburn_ZN/ZN_ssWSB.f90)
set(PY_SRCS
    Rocburn_PY/data_py.f90
    Rocburn_PY/init_py.f90
    Rocburn_PY/Rocburn_py.f90
    Rocburn_PY/setup_py.f90
    Rocburn_PY/update_py.f90)
set(LIB_SRCS
    ${CORE_SRCS}
    ${APN_SRCS}
    ${ZN_SRCS}
    ${PY_SRCS})

add_library(Rocburn ${LIB_SRCS})
target_link_libraries(Rocburn IMPACT::SITCOMF)
set_target_properties(Rocburn PROPERTIES VERSION ${ROCSTAR_VERSION}
        SOVERSION ${ROCSTAR_MAJOR_VERSION})

install(TARGETS Rocburn
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
