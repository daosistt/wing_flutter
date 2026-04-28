#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <gmsh.h>

#include "CalcMesh.h"

int main() {
    const unsigned int GMSH_TRIANGLE_CODE = 2;  // вместо GMSH_TETR_CODE = 4

    gmsh::initialize();
    gmsh::model::add("ellipse_model");

    try {
        gmsh::merge("../mesh/ellipse.msh");
    } catch(...) {
        gmsh::logger::write("Could not load MSH mesh: bye!");
        gmsh::finalize();
        return -1;
    }

    // Извлечём из gmsh данные об узлах сетки
    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    // И данные об элементах сетки - ищем треугольники
    std::vector<std::size_t>* trianglesNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] != GMSH_TRIANGLE_CODE)
            continue;
        trianglesNodesTags = &elementNodeTags[i];
    }

    if(trianglesNodesTags == nullptr) {
        std::cout << "Can not find triangle data. Exiting." << std::endl;
        gmsh::finalize();
        return -2;
    }

    std::cout << "The model has " << nodeTags.size() << " nodes and " 
              << trianglesNodesTags->size() / 3 << " triangles." << std::endl;

    // Проверим, что номера узлов идут подряд
    for(unsigned int i = 0; i < nodeTags.size(); ++i) {
        assert(i == nodeTags[i] - 1);
    }
    assert(trianglesNodesTags->size() % 3 == 0);  // 3 вершины на треугольник

    CalcMesh mesh(nodesCoord, *trianglesNodesTags);

    // Параметры эллипса:
    // a = 2.0 - большая полуось (по x)
    // b = 1.0 - малая полуось (по y)
    // Вращение вокруг малой полуоси (ось y)
    
    double angular_velocity = 2 * M_PI;  // полный оборот за 1 секунду

    mesh.setVelocities(
        // v_x: вращение вокруг y
        [angular_velocity](double x, double y, double z, double t) {
            return -angular_velocity * z;
        },
        // v_y: не меняется при вращении вокруг y
        [](double x, double y, double z, double t) {
            return 0.0;
        },
        // v_z: вращение вокруг y
        [angular_velocity](double x, double y, double z, double t) {
            return angular_velocity * x;
        }
    );

    mesh.setAngularVelocity(
        [](double t) { return 0.0; },
        [angular_velocity](double t) { return angular_velocity; },
        [](double t) { return 0.0; }
    );

    gmsh::finalize();

    mesh.snapshot(0);
    double tau = 0.01;
    for(unsigned int step = 1; step < 151; step++) {
        mesh.doTimeStep(tau);
        mesh.snapshot(step);
    }

    return 0;
}
