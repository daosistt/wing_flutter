#ifndef CALC_NODE_H
#define CALC_NODE_H

/**
 * class CalcNode - Узел расчётной сетки
 */
class CalcNode {
    friend class CalcMesh;

  protected:
    double x, y, z;             /**< Текущие абсолютные координаты узла */
    double relX, relY, relZ;    /**< Координаты относительно центра масс */

    double smth;                /**< Скалярное поле (например, температура/плотность) */
    double vx, vy, vz;          /**< Текущие скорости узла */


  public:
    /**
     * CalcNode - Конструктор по умолчанию
     */
    CalcNode();

    /**
     * CalcNode - Конструктор
     * @x: Начальная X-координата
     * @y: Начальная Y-координата
     * @z: Начальная Z-координата
     * @smth: Значение скалярного поля
     * @vx: Начальная скорость по X
     * @vy: Начальная скорость по Y
     * @vz: Начальная скорость по Z
     */
    CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz);

    /**
     * move - Переместить узел согласно его скорости
     * @tau: Шаг по времени
     */
    void move(double tau);
};

#endif // CALC_NODE_H
