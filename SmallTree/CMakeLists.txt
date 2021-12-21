
add_library(SmallTree STATIC f90_tools.f90 impedance_init_sub.f90 impedance_sub.f90 root_imp.f90 special_functions.f90)

target_include_directories(SmallTree
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        )