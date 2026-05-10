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

    double bendQ;              /**< q - Обобщенная координата изгиба: вертикальное смещение свободного конца крыла вдоль оси Z */
    double bendQDot;           /**< q' - Обобщенная скорость изгиба: скорость изменения bendQ */
    double bendQDDot;          /**< q'' - Обобщенное ускорение изгиба: ускорение свободного конца крыла вдоль оси Z */
    double bendMass;           /**< M - Эффективная масса изгибной формы, входит в M*q'' + C*q' + K*q = Q */
    double bendDamping;        /**< C - Коэффициент механического демпфирования изгиба, входит в член C*q' */
    double bendStiffness;      /**< K - Жесткость изгиба крыла, входит в член K*q */
    double lastBendingForce;   /**< Q - Последняя вычисленная обобщенная аэродинамическая сила изгиба Q */
    double wingSpanLength;     /**< L_wing - Длина крыла вдоль оси Y от закрепления до свободного края */
    double wingRootY;          /**< Координата закрепленного корня крыла по оси Y */
    double wingTipY;           /**< Координата свободного конца крыла по оси Y */

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

    /**
     * writePvd - создать Pvd-файл для VTK
     * @filename: название файла
     * @snapshots: 
     * @tau: момент времени (в секундах __вероятно__)
     */
    void writePvd(const std::string& filename, unsigned int snapshots, double tau) const;

    /**
     * writePvd - создать Pvd-файл для VTK
     * @windVx: скорость ветра
     * @windVy: угол при инициалиазации
     * @airDensity: плотность воздуха
     * @dragCoefficient: 
     * @angularDamping: угловое демпфирование 
     */
    void configureWindFlutter(double windVx,
                              double windVy,
                              double windVz,
                              double initialAngle,
                              double airDensity,
                              double dragCoefficient,
                              double angularDamping);

    /**
     * @initialBend: начальный сдвиг
     * @initialBendVelocity: начальная скорость сдвига
     * @bendMassValue: 
     * @bendStiffnessValue: жесткость
     * @bendDampingValue: демпфирование
     */
    void configureWingBending(double initialBend,
                          double initialBendVelocity,
                          double bendMassValue,
                          double bendStiffnessValue,
                          double bendDampingValue);

      /**< Функция формы изгиба крыла вдоль оси Y */
      double bendingShape(double relY) const; 

    /**
     * updateCenterOfMass - Пересчитать положение центра масс
     */
    void updateCenterOfMass();
};

#endif // CALC_MESH_H
