#ifndef CALC_MESH_H
#define CALC_MESH_H

#include <vector>
#include <string>
#include "CalcNode.h"
#include "Element.h"
#include <functional>

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
 * @yaw_rate_func: Функция, возвращающая угловую скорость рыскания в момент времени
 * @pitch_rate_func: Функция, возвращающая угловую скорость тангажа в момент времени
 * @roll_rate_func: Функция, возвращающая угловую скорость крена в момент времени
 */
struct CenterOfMass {
    double centerX, centerY, centerZ;
    double vx, vy, vz;
    
    double yaw, pitch, roll;
   
    std::function<double(double)> yaw_rate_func;
    std::function<double(double)> pitch_rate_func;
    std::function<double(double)> roll_rate_func;

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

    /** Функции поля скоростей (зависят от положения и времени) */
    std::function<double(double,double,double,double)> vxFunc;
    std::function<double(double,double,double,double)> vyFunc;
    std::function<double(double,double,double,double)> vzFunc;

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

    /**
     * setVelocities - Задать функции поля скоростей
     * @fx: Функция, возвращающая vx(x,y,z,t)
     * @fy: Функция, возвращающая vy(x,y,z,t)
     * @fz: Функция, возвращающая vz(x,y,z,t)
     */
    void setVelocities(std::function<double(double,double,double,double)> fx,
                       std::function<double(double,double,double,double)> fy,
                       std::function<double(double,double,double,double)> fz);

    /**
     * updateCenterOfMass - Пересчитать положение центра масс
     */
    void updateCenterOfMass();

    /**
     * setAngularVelocity - Задать функции угловых скоростей
     * @yaw_rate: Функция, возвращающая угловую скорость рыскания (вокруг Z)
     * @pitch_rate: Функция, возвращающая угловую скорость тангажа (вокруг Y)
     * @roll_rate: Функция, возвращающая угловую скорость крена (вокруг X)
     */
    void setAngularVelocity(
        std::function<double(double)> yaw_rate,
        std::function<double(double)> pitch_rate,
        std::function<double(double)> roll_rate
    );
};

#endif // CALC_MESH_H
