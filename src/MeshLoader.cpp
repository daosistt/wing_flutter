#include "MeshLoader.h"
#include <gmsh.h>
#include <iostream>
#include <cmath>

MeshLoader::MeshLoader(const std::string& filename) : filename(filename) {}

MeshLoader& MeshLoader::setParams(double size, int algo, bool optimize, double angle) {
    this->meshSize = size;
    this->algorithm = algo;
    this->optimize = optimize;
    this->angle = angle;
    return *this;
}

MeshData MeshLoader::generate() {
    const unsigned int GMSH_TETR_CODE = 4;

    gmsh::initialize();
    gmsh::model::add("mesh");

    try {
        gmsh::merge(filename);
    } catch(...) {
        gmsh::finalize();
        throw std::runtime_error("Could not load STL mesh: " + filename);
    }

    // Восстановление геометрии
    gmsh::option::setNumber("Mesh.StlRemoveDuplicateTriangles", 1);
    double angleRad = angle * M_PI / 180.0;
    gmsh::model::mesh::classifySurfaces(angleRad, true, true, angleRad);
    gmsh::model::mesh::createGeometry();

    // Создание объёма
    std::vector<std::pair<int, int>> surfaces;
    gmsh::model::getEntities(surfaces, 2);
    std::vector<int> surfaceTags;
    for (auto& s : surfaces) surfaceTags.push_back(s.second);
    int loop = gmsh::model::geo::addSurfaceLoop(surfaceTags);
    gmsh::model::geo::addVolume({loop});
    gmsh::model::geo::synchronize();

    // Настройки сетки
    gmsh::option::setNumber("Mesh.Algorithm3D", 4);
    gmsh::option::setNumber("Mesh.Algorithm", algorithm);
    if (optimize) {
        gmsh::option::setNumber("Mesh.OptimizeNetgen", 1);
        gmsh::option::setNumber("Mesh.Optimize", 1);
    }
    gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", 1);
    gmsh::option::setNumber("Mesh.MinimumCirclePoints", 20);

    // Поле размера
    int field = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(field, "F", std::to_string(meshSize));
    gmsh::model::mesh::field::setAsBackgroundMesh(field);

    // Генерация
    gmsh::model::mesh::generate(3);

    // Извлечение данных
    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    std::vector<std::size_t>* tetras = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);

    for (size_t i = 0; i < elementTypes.size(); ++i) {
        if (elementTypes[i] == GMSH_TETR_CODE) {
            tetras = &elementNodeTags[i];
            break;
        }
    }

    checkTetraData(tetras);

    // Конвертация в 0-индексацию
    std::vector<std::size_t> tetraNodes(tetras->size());
    for (size_t i = 0; i < tetras->size(); ++i) {
        tetraNodes[i] = (*tetras)[i] - 1;
    }

    gmsh::finalize();

    return {nodesCoord, tetraNodes};
}

void MeshLoader::checkTetraData(const std::vector<std::size_t>* tetras) const {
    if (!tetras) {
        throw std::runtime_error("No tetrahedra found in mesh");
    }
    if (tetras->size() % 4 != 0) {
        throw std::runtime_error("Invalid tetra data size");
    }
}
