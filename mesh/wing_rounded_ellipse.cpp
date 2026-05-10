#include <gmsh.h>

#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <algorithm>

struct Node {
    double x, y, z;
};

int main(int argc, char **argv)
{
    gmsh::initialize();

    gmsh::model::add("rounded_ellipse_wing");

    // -----------------------------
    // Геометрия крыла
    // -----------------------------
    // Оси:
    // x — ширина крыла / хорда, вдоль ветра
    // y — длина крыла, от корня к свободному концу
    // z — высота / толщина

    double chord = 1.0;   // ширина крыла по x
    double span  = 4.0;   // длина крыла по y
    double T     = 0.12;  // максимальная толщина по z
    double R     = 0.25;  // радиус скругления концов в плане x-y

    // -----------------------------
    // Параметры сетки
    // -----------------------------
    int Ny = 40;          // по длине y
    int Nx = 12;          // по ширине x
    int Nz = 1;           // по толщине z: толщина в 1 элемент

    // Чтобы не получить нулевую толщину на самой кромке эллипса
    double etaMax = 0.97;

    // Формат msh
    gmsh::option::setNumber("Mesh.MshFileVersion", 4.1);
    // Для старого формата можно заменить на:
    // gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);

    std::vector<Node> nodes;
    std::vector<std::array<int, 4>> tets;

    auto halfChord = [&](double y) {
        double h = chord / 2.0;

        // Скругление у корня y = 0
        if (y < R) {
            double dy = y - R;
            h = chord / 2.0 - R
              + std::sqrt(std::max(0.0, R * R - dy * dy));
        }

        // Скругление у свободного конца y = span
        else if (y > span - R) {
            double dy = y - (span - R);
            h = chord / 2.0 - R
              + std::sqrt(std::max(0.0, R * R - dy * dy));
        }

        return h;
    };

    auto nodeId = [&](int iy, int ix, int iz) {
        return iz + (Nz + 1) * (ix + (Nx + 1) * iy);
    };

    // -----------------------------
    // Генерация узлов
    // -----------------------------
    for (int iy = 0; iy <= Ny; ++iy) {

        double y = span * double(iy) / double(Ny);
        double halfX = halfChord(y);

        for (int ix = 0; ix <= Nx; ++ix) {

            double eta = -etaMax
                       + 2.0 * etaMax * double(ix) / double(Nx);

            double x = eta * halfX;

            // Эллиптическое поперечное сечение x-z
            double zMax = 0.5 * T
                        * std::sqrt(std::max(0.0, 1.0 - eta * eta));

            for (int iz = 0; iz <= Nz; ++iz) {

                double s = -1.0
                         + 2.0 * double(iz) / double(Nz);

                double z = s * zMax;

                nodes.push_back({x, y, z});
            }
        }
    }

    // -----------------------------
    // Генерация тетраэдров
    // -----------------------------
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            for (int iz = 0; iz < Nz; ++iz) {

                int n000 = nodeId(iy,     ix,     iz);
                int n100 = nodeId(iy + 1, ix,     iz);
                int n010 = nodeId(iy,     ix + 1, iz);
                int n110 = nodeId(iy + 1, ix + 1, iz);

                int n001 = nodeId(iy,     ix,     iz + 1);
                int n101 = nodeId(iy + 1, ix,     iz + 1);
                int n011 = nodeId(iy,     ix + 1, iz + 1);
                int n111 = nodeId(iy + 1, ix + 1, iz + 1);

                tets.push_back({n000, n100, n010, n001});
                tets.push_back({n100, n110, n010, n111});
                tets.push_back({n100, n010, n001, n111});
                tets.push_back({n100, n001, n101, n111});
                tets.push_back({n010, n001, n011, n111});
            }
        }
    }

    // -----------------------------
    // Передача сетки в Gmsh
    // -----------------------------
    int volumeTag = 1;

    gmsh::model::addDiscreteEntity(3, volumeTag);

    std::vector<std::size_t> nodeTags;
    std::vector<double> coords;

    nodeTags.reserve(nodes.size());
    coords.reserve(3 * nodes.size());

    for (std::size_t i = 0; i < nodes.size(); ++i) {
        nodeTags.push_back(i + 1);

        coords.push_back(nodes[i].x);
        coords.push_back(nodes[i].y);
        coords.push_back(nodes[i].z);
    }

    gmsh::model::mesh::addNodes(
        3,
        volumeTag,
        nodeTags,
        coords
    );

    std::vector<std::size_t> elemTags;
    std::vector<std::size_t> elemNodeTags;

    elemTags.reserve(tets.size());
    elemNodeTags.reserve(4 * tets.size());

    for (std::size_t e = 0; e < tets.size(); ++e) {
        elemTags.push_back(e + 1);

        elemNodeTags.push_back(tets[e][0] + 1);
        elemNodeTags.push_back(tets[e][1] + 1);
        elemNodeTags.push_back(tets[e][2] + 1);
        elemNodeTags.push_back(tets[e][3] + 1);
    }

    // Gmsh element type:
    // 4 = 4-node tetrahedron
    gmsh::model::mesh::addElementsByType(
        volumeTag,
        4,
        elemTags,
        elemNodeTags
    );

    // -----------------------------
    // Физические группы
    // -----------------------------
    gmsh::model::addPhysicalGroup(3, {volumeTag}, 1);
    gmsh::model::setPhysicalName(3, 1, "wing_volume");

    // -----------------------------
    // Запись msh
    // -----------------------------
    gmsh::write("wing_rounded_ellipse.msh");

    gmsh::finalize();

    std::cout << "Written: wing_rounded_ellipse.msh\n";
    std::cout << "Nodes: " << nodes.size() << "\n";
    std::cout << "Tetrahedra: " << tets.size() << "\n";

    return 0;
}