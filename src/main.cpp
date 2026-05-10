#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <gmsh.h>

#include "CalcMesh.h"

int main() {
    const unsigned int GMSH_TETR_CODE = 4;

    gmsh::initialize();
    gmsh::model::add("ellipse_model");

    try {
        // gmsh::merge("../mesh/ellipse.msh");
        gmsh::merge("../mesh/wing_rounded_ellipse.msh");
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

    // И данные об элементах сетки - ищем тетраэдры
    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] != GMSH_TETR_CODE)
            continue;
        tetrsNodesTags = &elementNodeTags[i];
    }

    if(tetrsNodesTags == nullptr) {
        std::cout << "Can not find tetra data. Exiting." << std::endl;
        gmsh::finalize();
        return -2;
    }

    std::cout << "The model has " << nodeTags.size() << " nodes and " 
              << tetrsNodesTags->size() / 4 << " tetrs." << std::endl;

    // Проверим, что номера узлов идут подряд
    for(unsigned int i = 0; i < nodeTags.size(); ++i) {
        assert(i == nodeTags[i] - 1);
    }
    assert(tetrsNodesTags->size() % 4 == 0);

    CalcMesh mesh(nodesCoord, *tetrsNodesTags);

    gmsh::finalize();

    // Reduced-order FSI model:
    // The wind is directed along +Y, perpendicular to the initial major axis.
    // The body is hinged at the rightmost point of the major semi-axis and
    // rotates in the XY plane around the Z axis.
    const double windVx = 8.0;
    const double windVy = 0.0;
    const double windVz = 8.0;
    const double initialAngle = 0.0;
    const double airDensity = 1.225;
    const double dragCoefficient = 1.15;
    const double angularDamping = 0.04;
    mesh.configureWindFlutter(windVx, windVy, windVz, initialAngle,
                              airDensity, dragCoefficient, angularDamping);

    const double initialBend = 0.0;
    const double initialBendVelocity = 0.0;

    const double bendMass = 0.25;
    const double bendStiffness = 20.0;
    const double bendDamping = 0.2;

    mesh.configureWingBending(initialBend,
                              initialBendVelocity,
                              bendMass,
                              bendStiffness,
                              bendDamping);

    mesh.snapshot(0);
    const double tau = 0.01;
    const unsigned int steps = 300;
    for(unsigned int step = 1; step <= steps; step++) {
        mesh.doTimeStep(tau);
        mesh.snapshot(step);
    }
    mesh.writePvd("wing.pvd", steps + 1, tau);

    return 0;
}
