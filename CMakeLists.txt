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

# VTK
find_package(VTK REQUIRED)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${VTK_INCLUDE_DIRS}
)

add_library(mesh_lib
    src/CalcNode.cpp
    src/CalcMesh.cpp
    src/MeshLoader.cpp
)

add_executable(mv_wing src/main.cpp)

target_link_libraries(mesh_lib
    ${VTK_LIBRARIES}
    gmsh
)

target_link_libraries(mv_wing
    mesh_lib
)