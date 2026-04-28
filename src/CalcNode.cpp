#include "CalcNode.h"

CalcNode::CalcNode() : x(0.0), y(0.0), z(0.0), smth(0.0), vx(0.0), vy(0.0), vz(0.0)
{
}

CalcNode::CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz) 
        : x(x), y(y), z(z), smth(smth), vx(vx), vy(vy), vz(vz), relX(x), relY(y), relZ(z)
{
}

void CalcNode::move(double tau) {
    x += vx * tau;
    y += vy * tau;
    z += vz * tau;
}
