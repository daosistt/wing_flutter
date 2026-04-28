#ifndef ELEMENT_H
#define ELEMENT_H

/**
 * class Element - Тетраэдральный элемент сетки
 */
class Element {
    friend class CalcMesh;

  protected:
    unsigned long nodesIds[4];  /**< Индексы узлов, образующих тетраэдр */
};

#endif // ELEMENT_H
