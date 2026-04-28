#ifndef MESH_LOADER_H
#define MESH_LOADER_H

#include <string>
#include <vector>
#include <tuple>

/**
 * struct MeshData - Результат загрузки и построения сетки
 * @nodes: Вектор координат узлов [x0,y0,z0, x1,y1,z1, ...]
 * @tetraNodes: Вектор индексов узлов тетраэдров (0-индексация)
 */
struct MeshData {
    std::vector<double> nodes;
    std::vector<std::size_t> tetraNodes;
};



/**
 * class MeshLoader - Загрузчик STL и построитель тетраэдральной сетки
 */
class MeshLoader {
public:
    /**
     * MeshLoader - Конструктор
     * @filename: Путь к STL-файлу
     */
    MeshLoader(const std::string& filename);


    /**
     * setParams - Установить параметры построения сетки
     * @size: Характерный размер элемента
     * @algo: Алгоритм построения сетки (по умолчанию 6)
     * @optimize: Включить оптимизацию сетки (true/false)
     * @angle: Угол для классификации поверхностей (градусы)
     */
    MeshLoader& setParams(double size, int algo, bool optimize, double angle);

    /**
     * generate - Построить сетку и извлечь данные
     *
     * @return MeshData с узлами и тетраэдрами
     * @throws std::runtime_error при ошибке
     */
    MeshData generate();

private:
    std::string filename;
    double meshSize = 16.0;
    int algorithm = 6;
    bool optimize = true;
    double angle = 40.0;

    void checkTetraData(const std::vector<std::size_t>* tetras) const;
};

#endif // MESH_LOADER_h
