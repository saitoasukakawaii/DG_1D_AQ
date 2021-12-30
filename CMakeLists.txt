cmake_minimum_required(VERSION 3.16)
project(DG_1D_AQ)


enable_language(Fortran)
set(INC_DIR ./ ~/eigen-3.4.0 ~/boost_1_77_0 /usr/local/include)

set(CMAKE_CXX_FLAGS "-Wall -std=c++17")
set(CMAKE_CXX_FLAGS_DEBUG "-D_DEBUG -O0 -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

include_directories(${INC_DIR})

add_subdirectory(Polynomial)
list(APPEND EXTRA_LIBS Polynomial)

add_subdirectory(SmallTree)
list(APPEND EXTRA_LIBS SmallTree)

list(APPEND EXTRA_LIBS lapack blas fftw3)

add_executable(DG_1D_AQ main.cpp Element.cpp artery.cpp tools.h)

target_link_libraries(DG_1D_AQ PUBLIC ${EXTRA_LIBS})

target_include_directories(DG_1D_AQ PUBLIC
        "${PROJECT_BINARY_DIR}"
        )

add_custom_command(
        TARGET DG_1D_AQ POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_SOURCE_DIR}/input.dat
        ${CMAKE_CURRENT_BINARY_DIR}/input.dat
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_SOURCE_DIR}/Geometry/topology55.txt
        ${CMAKE_CURRENT_BINARY_DIR}/Geometry/topology55.txt
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_SOURCE_DIR}/Geometry/geometry55.txt
        ${CMAKE_CURRENT_BINARY_DIR}/Geometry/geometry55.txt)
