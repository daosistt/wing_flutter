cmake_minimum_required(VERSION 3.10)
project(WingMesh)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(VTK_LIBRARIES
    vtkCommonCore
    vtkCommonDataModel
    vtkIOCore
    vtkIOXML
    vtksys
)

set(GMSH_LIBRARIES gmsh)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${VTK_INCLUDE_DIRS}
    ${GMSH_INCLUDE_DIRS}
)

add_library(mesh_lib
    src/CalcNode.cpp
    src/CalcMesh.cpp
    src/MeshLoader.cpp
)

# Создаем исполняемый файл
add_executable(mv_wing src/main.cpp)

# Линкуем библиотеки
target_link_libraries(mesh_lib ${VTK_LIBRARIES})
target_link_libraries(mv_wing mesh_lib ${GMSH_LIBRARIES})
