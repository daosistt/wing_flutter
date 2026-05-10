cmake_minimum_required(VERSION 3.10)
project(WingMesh)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# GMSH
add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/opt/homebrew/lib/libgmsh.dylib"
    INTERFACE_INCLUDE_DIRECTORIES "/opt/homebrew/include"
)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_executable(wing_rounded_ellipse wing_rounded_ellipse.cpp)

target_link_libraries(wing_rounded_ellipse
    gmsh
)