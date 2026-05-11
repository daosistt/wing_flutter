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
    // ВЕТЕР
    const double windVx = 20.0;
    const double windVy = 0.0;
    const double windVz = 0.0;
    const double initialAngle = 0.0;
    const double airDensity = 1.225;
    const double dragCoefficient = 1.15;
    const double angularDamping = 0.04;
    mesh.configureWindFlutter(windVx, windVy, windVz, initialAngle,
                              airDensity, dragCoefficient, angularDamping);

    const double liftSlope = 10;
    const double maxAlphaEff = 0.35;
    mesh.configureFlutterAerodynamics(liftSlope, maxAlphaEff);

    const double alphataulag = 0.015;
    mesh.configureAeroLag(false, alphataulag);

    // СДВИГ
    const double initialBend = 0.1;
    const double initialBendVelocity = 0.0;

    const double bendMass = 0.25;
    const double bendStiffness = 40.0;
    const double bendDamping = 0.001;

    mesh.configureWingBending(initialBend,
                              initialBendVelocity,
                              bendMass,
                              bendStiffness,
                              bendDamping);
    // КРУЧЕНИЕ
    const double initialTorsion = 0.01;
    const double initialTorsionVelocity = 0.0;

    const double torsionInertia = 0.01;
    const double torsionStiffness = 100;
    const double torsionDamping = 0.9;

    mesh.configureWingTorsion(initialTorsion,
                            initialTorsionVelocity,
                            torsionInertia,
                            torsionStiffness,
                            torsionDamping);

    // МАТЕРИАЛ
    const double youngModulus = 70.0e9;    // Па, алюминий
    const double poissonRatio = 0.33;      // безразмерный
    const double yieldStress = 300.0e6;    // Па

    mesh.configureWingMaterial(youngModulus,
                            poissonRatio,
                            yieldStress);


    mesh.snapshot(0);
    const double tau = 0.002;

    mesh.configureTau(tau);
    const unsigned int steps = 3000;
    for(unsigned int step = 1; step <= steps; step++) {
        mesh.doTimeStep();
        mesh.snapshot(step);
    }
    mesh.writePvd("wing.pvd", steps + 1, tau);

    return 0;
}
