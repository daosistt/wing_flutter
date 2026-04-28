#include <set>
#include <string>
#include <gmsh.h>
#include <cmath>

int main(int argc, char **argv)
{
    gmsh::initialize();

    gmsh::model::add("ellipse");

    double lc = 1e-1;
    double a  = 2.0;   // большая полуось (по x)
    double b  = 1.0;   // малая полуось (по y)

    // Центр
    gmsh::model::geo::addPoint(0.0,  0.0, 0.0, lc, 1);

    // Точки на эллипсе (для дуг) - вершины эллипса
    gmsh::model::geo::addPoint( a,   0.0, 0.0, lc, 2);  // правая вершина
    gmsh::model::geo::addPoint( 0.0,  b,  0.0, lc, 3);  // верхняя вершина
    gmsh::model::geo::addPoint(-a,   0.0, 0.0, lc, 4);  // левая вершина
    gmsh::model::geo::addPoint( 0.0, -b,  0.0, lc, 5);  // нижняя вершина

    gmsh::model::geo::addEllipseArc(2, 1, 4, 3, 1);  // от правой вершины к верхней
    gmsh::model::geo::addEllipseArc(3, 1, 4, 4, 2);  // от верхней вершины к левой
    gmsh::model::geo::addEllipseArc(4, 1, 2, 5, 3);  // от левой вершины к нижней
    gmsh::model::geo::addEllipseArc(5, 1, 2, 2, 4);  // от нижней вершины к правой

    // Контур и поверхность
    gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);

    gmsh::model::geo::synchronize();

    // Граница (для BC Дирихле)
    gmsh::model::addPhysicalGroup(1, {1, 2, 3, 4}, 1);
    gmsh::model::setPhysicalName(1, 1, "boundary");

    // Объём (поверхность)
    gmsh::model::addPhysicalGroup(2, {1}, 2);
    gmsh::model::setPhysicalName(2, 2, "domain");

    gmsh::model::mesh::generate(2);

    gmsh::write("ellipse.msh");

    gmsh::finalize();
    return 0;
}
