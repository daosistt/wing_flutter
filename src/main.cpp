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

    // Параметры эллипса:
    // a = 2.0 - большая полуось (по x)
    // b = 1.0 - малая полуось (по y)
    // Колебания: вращение вокруг оси y с переменной угловой скоростью
    
    double amplitude = 0.2 * M_PI;      // амплитуда угловой скорости
    double frequency = 0.5;           // частота колебаний (Гц)
    
    // Функция угловой скорости от времени (синусоидальные колебания)
    auto omega = [amplitude, frequency](double t) -> double {
        return amplitude * sin(2 * M_PI * frequency * t);
    };

    mesh.setVelocities(
        // v_x: зависит от угловой скорости в момент t
        [omega](double x, double y, double z, double t) {
            return -omega(t) * z;
        },
        // v_y: не меняется
        [](double x, double y, double z, double t) {
            return 0.0;
        },
        // v_z: зависит от угловой скорости в момент t
        [omega](double x, double y, double z, double t) {
            return omega(t) * x;
        }
    );

    mesh.setAngularVelocity(
        [](double t) { return 0.0; },
        omega,
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
