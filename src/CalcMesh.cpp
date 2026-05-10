#include "CalcMesh.h"
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>

namespace {
const double PI = 3.14159265358979323846;

double sqr(double value) {
    return value * value;
}
}

CalcMesh::CalcMesh(const std::vector<double>& nodesCoords,
                   const std::vector<std::size_t>& tetrsPoints)
    : currentTime(0.0),
      hingeX(0.0), hingeY(0.0), hingeZ(0.0),

      angle(0.0), angularVelocity(0.0), lastTorque(0.0),
      inertiaZ(1.0), bodyMass(1.0), nodeArea(1.0),

      windVx(0.0), windVy(8.0), windVz(0.0),
      airDensity(1.225), dragCoefficient(1.15),
      angularDamping(0.35), referenceArea(0.0),

      bendQ(0.0), bendQDot(0.0), bendQDDot(0.0),
      bendMass(1.0), bendDamping(0.1), bendStiffness(10.0),
      lastBendingForce(0.0),
      wingSpanLength(1.0),
      wingRootY(0.0),
      wingTipY(0.0)
    {
    nodes.resize(nodesCoords.size() / 3);
    for(unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
        double pointX = nodesCoords[i*3];
        double pointY = nodesCoords[i*3 + 1];
        double pointZ = nodesCoords[i*3 + 2];

        double smth = 0.0;
        nodes[i] = CalcNode(pointX, pointY, pointZ, smth, 0.0, 0.0, 0.0);
    }

    double minX = nodes[0].x, maxX = nodes[0].x;
    double minY = nodes[0].y, maxY = nodes[0].y;
    double minZ = nodes[0].z, maxZ = nodes[0].z;
    for (const auto& node : nodes) {
        minX = std::min(minX, node.x);
        maxX = std::max(maxX, node.x);
        minY = std::min(minY, node.y);
        maxY = std::max(maxY, node.y);
        minZ = std::min(minZ, node.z);
        maxZ = std::max(maxZ, node.z);
    }

    hingeX = 0.5 * (minX + maxX);
    hingeY = maxY;
    hingeZ = 0.5 * (minZ + maxZ);

    wingRootY = hingeY;
    wingTipY = minY;
    wingSpanLength = std::max(wingRootY - wingTipY, 1.0e-9);

    for (auto& node : nodes) {
        node.relX = node.x - hingeX;
        node.relY = node.y - hingeY;
        node.relZ = node.z - hingeZ;
    }

    const double a = 0.5 * (maxX - minX);
    const double b = 0.5 * (maxY - minY);
    const double thickness = std::max(maxZ - minZ, 1.0e-6);
    bodyMass = PI * a * b * thickness;
    referenceArea = std::max(2.0 * b * thickness, 1.0e-6);
    nodeArea = referenceArea / static_cast<double>(nodes.size());

    inertiaZ = 0.0;
    const double pointMass = bodyMass / static_cast<double>(nodes.size());
    for (const auto& node : nodes) {
        inertiaZ += pointMass * (sqr(node.relX) + sqr(node.relY));
    }
    inertiaZ = std::max(inertiaZ, 1.0e-6);

    elements.resize(tetrsPoints.size() / 4);
    for(unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
        elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
        elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
        elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
        elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
    }

    updateKinematics();
    updateAerodynamicFields();
}

void CalcMesh::doTimeStep(double tau) {
    /*
     * 1. Считаем аэродинамические силы при текущем положении
     *    и текущих скоростях узлов.
     */
    updateAerodynamicFields();

    /*
     * 2. Уравнение одной изгибной степени свободы:
     *
     * bendMass * bendQDDot
     * + bendDamping * bendQDot
     * + bendStiffness * bendQ
     * = lastBendingForce
     */
    bendQDDot =
        (lastBendingForce
         - bendDamping * bendQDot
         - bendStiffness * bendQ)
        / bendMass;

    /*
     * 3. Semi-implicit Euler:
     *
     * Сначала обновляем скорость,
     * потом положение.
     *
     * Это устойчивее обычного явного Эйлера.
     */
    bendQDot += bendQDDot * tau;
    bendQ    += bendQDot  * tau;

    /*
     * 4. Обновляем координаты узлов по новой амплитуде изгиба.
     */
    updateKinematics();

    /*
     * 5. Пересчитываем силы для записи в VTK на новом положении.
     */
    updateAerodynamicFields();

    currentTime += tau;
}

void CalcMesh::configureWindFlutter(double windVxValue,
                                    double windVyValue,
                                    double windVzValue,
                                    double initialAngle,
                                    double airDensityValue,
                                    double dragCoefficientValue,
                                    double angularDampingValue) {
    windVx = windVxValue;
    windVy = windVyValue;
    windVz = windVzValue;
    // windVz = 0.0;
    angle = initialAngle;
    airDensity = airDensityValue;
    dragCoefficient = dragCoefficientValue;
    angularDamping = angularDampingValue;
    angularVelocity = 0.0;
    currentTime = 0.0;
    updateKinematics();
    updateAerodynamicFields();
}

void CalcMesh::configureWingBending(double initialBend,
                                    double initialBendVelocity,
                                    double bendMassValue,
                                    double bendStiffnessValue,
                                    double bendDampingValue) {
    bendQ = initialBend;
    bendQDot = initialBendVelocity;
    bendQDDot = 0.0;

    bendMass = std::max(bendMassValue, 1.0e-9);
    bendStiffness = std::max(bendStiffnessValue, 0.0);
    bendDamping = std::max(bendDampingValue, 0.0);

    lastBendingForce = 0.0;

    updateKinematics();
    updateAerodynamicFields();
}

double CalcMesh::bendingShape(double relY) const {
    /*
     * relY = node.y0 - hingeY.
     *
     * Корень крыла:
     *   relY = 0
     *   xi = 0
     *
     * Свободный конец:
     *   relY = -wingSpanLength
     *   xi = 1
     */

    double xi = -relY / wingSpanLength;

    xi = std::clamp(xi, 0.0, 1.0);

    /*
     * Консольная форма:
     * phi(0) = 0  - у закрепления нет смещения
     * phi'(0)= 0  - у закрепления нет поворота
     * phi(1) = 1  - bendQ равен смещению свободного конца
     */
    return xi * xi * (3.0 - 2.0 * xi);
}

void CalcMesh::updateKinematics() {
    com.centerX = 0.0;
    com.centerY = 0.0;
    com.centerZ = 0.0;

    for (auto& node : nodes) {
        const double phi = bendingShape(node.relY);

        /*
         * X не меняем:
         * вдоль X дует ветер.
         */
        node.x = hingeX + node.relX;

        /*
         * Y не меняем:
         * вдоль Y расположен размах крыла.
         */
        node.y = hingeY + node.relY;

        /*
         * Z меняем:
         * это вертикальный изгиб крыла.
         */
        node.z = hingeZ + node.relZ + bendQ * phi;

        /*
         * Скорость узла появляется только по Z,
         * потому что bendQ зависит от времени.
         */
        node.vx = 0.0;
        node.vy = 0.0;
        node.vz = bendQDot * phi;

        com.centerX += node.x;
        com.centerY += node.y;
        com.centerZ += node.z;
    }

    const double invN = 1.0 / static_cast<double>(nodes.size());

    com.centerX *= invN;
    com.centerY *= invN;
    com.centerZ *= invN;
}

void CalcMesh::updateAerodynamicFields() {
    lastTorque = 0.0;
    lastBendingForce = 0.0;

    const double windNorm = std::sqrt(sqr(windVx) + sqr(windVy) + sqr(windVz));
    const double windDirX = windNorm > 1.0e-9 ? windVx / windNorm : 0.0;
    const double windDirY = windNorm > 1.0e-9 ? windVy / windNorm : 0.0;
    double minProjection = nodes[0].x * windDirX + nodes[0].y * windDirY;
    double maxProjection = minProjection;
    for (const auto& node : nodes) {
        const double projection = node.x * windDirX + node.y * windDirY;
        minProjection = std::min(minProjection, projection);
        maxProjection = std::max(maxProjection, projection);
    }
    const double invSpan = 1.0 / std::max(maxProjection - minProjection, 1.0e-9);

    for(auto& node : nodes) {
        const double relFlowX = windVx - node.vx;
        const double relFlowY = windVy - node.vy;
        const double relFlowZ = windVz - node.vz;
        const double relSpeed = std::sqrt(sqr(relFlowX) + sqr(relFlowY) + sqr(relFlowZ));
        const double q = 0.5 * airDensity * relSpeed * relSpeed;
        const double projection = node.x * windDirX + node.y * windDirY;
        const double windward = (maxProjection - projection) * invSpan;
        const double pressureShape = 0.15 + 1.35 * windward * windward;
        const double localPressure = q * pressureShape;
        const double invSpeed = relSpeed > 1.0e-9 ? 1.0 / relSpeed : 0.0;

        node.pressure = localPressure;
        node.fx = dragCoefficient * localPressure * nodeArea * relFlowX * invSpeed;
        node.fy = dragCoefficient * localPressure * nodeArea * relFlowY * invSpeed;
        node.fz = dragCoefficient * localPressure * nodeArea * relFlowZ * invSpeed;

        const double phi = bendingShape(node.relY);

        /*
        * Обобщенная сила изгиба:
        *
        * z_i = z0_i + bendQ * phi_i
        *
        * dz_i / d(bendQ) = phi_i
        *
        * Поэтому:
        * Q = sum(Fz_i * phi_i)
        */
        lastBendingForce += node.fz * phi;

        const double rx = node.x - hingeX;
        const double ry = node.y - hingeY;
        lastTorque += rx * node.fy - ry * node.fx;
        node.smth = node.pressure;
    }
}



void CalcMesh::snapshot(unsigned int snap_number) {
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    auto pressure = vtkSmartPointer<vtkDoubleArray>::New();
    pressure->SetName("pressure");

    auto vel = vtkSmartPointer<vtkDoubleArray>::New();
    vel->SetName("velocity");
    vel->SetNumberOfComponents(3);

    auto force = vtkSmartPointer<vtkDoubleArray>::New();
    force->SetName("force");
    force->SetNumberOfComponents(3);

    auto wind = vtkSmartPointer<vtkDoubleArray>::New();
    wind->SetName("wind_velocity");
    wind->SetNumberOfComponents(3);

    auto displacement = vtkSmartPointer<vtkDoubleArray>::New();
    displacement->SetName("hinge_displacement");
    displacement->SetNumberOfComponents(3);

    for(unsigned int i = 0; i < nodes.size(); i++) {
        dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);
        double _vel[3] = {nodes[i].vx, nodes[i].vy, nodes[i].vz};
        double _force[3] = {nodes[i].fx, nodes[i].fy, nodes[i].fz};
        double _wind[3] = {windVx, windVy, windVz};
        double _displacement[3] = {
            nodes[i].x - (hingeX + nodes[i].relX),
            nodes[i].y - (hingeY + nodes[i].relY),
            0.0
        };
        vel->InsertNextTuple(_vel);
        force->InsertNextTuple(_force);
        wind->InsertNextTuple(_wind);
        displacement->InsertNextTuple(_displacement);
        pressure->InsertNextValue(nodes[i].pressure);
    }

    unstructuredGrid->SetPoints(dumpPoints);
    unstructuredGrid->GetPointData()->AddArray(vel);
    unstructuredGrid->GetPointData()->AddArray(force);
    unstructuredGrid->GetPointData()->AddArray(wind);
    unstructuredGrid->GetPointData()->AddArray(displacement);
    unstructuredGrid->GetPointData()->AddArray(pressure);

    auto timeArray = vtkSmartPointer<vtkDoubleArray>::New();
    timeArray->SetName("time");
    timeArray->InsertNextValue(currentTime);
    unstructuredGrid->GetFieldData()->AddArray(timeArray);

    auto stateArray = vtkSmartPointer<vtkDoubleArray>::New();
    stateArray->SetName("body_state");
    stateArray->SetNumberOfComponents(3);
    // double state[3] = {angle, angularVelocity, lastTorque};
    // stateArray->InsertNextTuple(state);
    unstructuredGrid->GetFieldData()->AddArray(stateArray);

    for(unsigned int i = 0; i < elements.size(); i++) {
        auto tetra = vtkSmartPointer<vtkTetra>::New();
        tetra->GetPointIds()->SetId(0, elements[i].nodesIds[0]);
        tetra->GetPointIds()->SetId(1, elements[i].nodesIds[1]);
        tetra->GetPointIds()->SetId(2, elements[i].nodesIds[2]);
        tetra->GetPointIds()->SetId(3, elements[i].nodesIds[3]);
        unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
    }

    std::string fileName = "wing-step-" + std::to_string(snap_number) + ".vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}

void CalcMesh::writePvd(const std::string& filename, unsigned int snapshots, double tau) const {
    std::ofstream out(filename.c_str());
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <Collection>\n";
    for (unsigned int step = 0; step < snapshots; ++step) {
        out << "    <DataSet timestep=\"" << step * tau
            << "\" group=\"\" part=\"0\" file=\"wing-step-" << step << ".vtu\"/>\n";
    }
    out << "  </Collection>\n";
    out << "</VTKFile>\n";
}

void CalcMesh::updateCenterOfMass() {
    com.centerX = com.centerY = com.centerZ = 0.0;
    for (auto& node : nodes) {
        com.centerX += node.x;
        com.centerY += node.y;
        com.centerZ += node.z;
    }
    com.centerX /= nodes.size();
    com.centerY /= nodes.size();
    com.centerZ /= nodes.size();
}
