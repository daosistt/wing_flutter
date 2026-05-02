#ifndef CALC_MESH_H
#define CALC_MESH_H

#include <vector>
#include <string>
#include "CalcNode.h"
#include "Element.h"

/**
 * struct CenterOfMass - Состояние центра масс сетки
 * @centerX: X-координата центра масс
 * @centerY: Y-координата центра масс
 * @centerZ: Z-координата центра масс
 * @vx: Линейная скорость по оси X
 * @vy: Линейная скорость по оси Y
 * @vz: Линейная скорость по оси Z
 * @yaw: Угол рыскания (вращение вокруг Z) в радианах
 * @pitch: Угол тангажа (вращение вокруг Y) в радианах
 * @roll: Угол крена (вращение вокруг X) в радианах
 */
struct CenterOfMass {
    double centerX, centerY, centerZ;
    double vx, vy, vz;
    
    double yaw, pitch, roll;

    CenterOfMass() 
        : centerX(0.0), centerY(0.0), centerZ(0.0)
        , vx(0.0), vy(0.0), vz(0.0)
        , yaw(0.0), pitch(0.0), roll(0.0) {}
};





/**
 * class CalcMesh - Основной класс сетки для моделирования тетраэдральных сеток
 *
 * Класс управляет тетраэдральной сеткой, выполняет интегрирование по времени,
 * применяет скорости и вращения, а также экспортирует состояния в формат VTK.
 */
class CalcMesh {
  protected:
    std::vector<CalcNode> nodes;    /**< Узлы сетки с их состоянием */
    std::vector<Element> elements;  /**< Тетраэдральные элементы */
    double currentTime;             /**< Текущее время моделирования */

    CenterOfMass com;               /**< Состояние центра масс */
    double hingeX, hingeY, hingeZ;   /**< Неподвижная точка подвеса */
    double angle;                    /**< Угол поворота в плоскости XY вокруг оси Z */
    double angularVelocity;          /**< Угловая скорость вокруг оси Z */
    double lastTorque;               /**< Момент аэродинамических сил вокруг подвеса */
    double inertiaZ;                 /**< Момент инерции относительно подвеса */
    double bodyMass;                 /**< Масса жёсткого тела в безразмерной модели */
    double nodeArea;                 /**< Эффективная площадь, отнесённая к одному узлу */
    double windVx, windVy, windVz;   /**< Скорость ветра */
    double airDensity;               /**< Плотность воздуха */
    double dragCoefficient;          /**< Коэффициент сопротивления */
    double angularDamping;           /**< Вязкое демпфирование вращения */
    double referenceArea;            /**< Характерная площадь обдува */

    void updateKinematics();
    void updateAerodynamicFields();

  public:
    /**
     * CalcMesh - Конструктор
     * @nodesCoords: Вектор координат узлов [x0,y0,z0, x1,y1,z1, ...]
     * @tetrsPoints: Вектор индексов узлов тетраэдров (индексация с 1 от Gmsh)
     */
    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints);

    /**
     * doTimeStep - Выполнить один шаг интегрирования по времени
     * @tau: Величина шага по времени (дельта t)
     */
    void doTimeStep(double tau);

    /**
     * snapshot - Сохранить текущее состояние сетки в VTK-файл
     * @snap_number: Номер состояния (используется в имени файла)
     */
    void snapshot(unsigned int snap_number);
    void writePvd(const std::string& filename, unsigned int snapshots, double tau) const;
    void configureWindFlutter(double windVx,
                              double windVy,
                              double initialAngle,
                              double airDensity,
                              double dragCoefficient,
                              double angularDamping);

    /**
     * updateCenterOfMass - Пересчитать положение центра масс
     */
    void updateCenterOfMass();
};

#endif // CALC_MESH_H
