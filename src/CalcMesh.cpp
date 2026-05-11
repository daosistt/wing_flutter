#include "CalcMesh.h"
#include "CalcNode.h"
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkCellData.h>
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
      wingTipY(0.0), 

      torsionQ(0.0), torsionQDot(0.0), torsionQDDot(0.0),
      torsionInertia(1.0), torsionDamping(0.05), torsionStiffness(5.0),
      lastTorsionMoment(0.0),
      torsionAxisX(0.0), torsionAxisZ(0.0),

      youngModulus(70.0e9), shearModulus(26.0e9),
      poissonRatio(0.33),
      yieldStress(300.0e6), maxEquivalentStress(0.0)
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

    torsionAxisX = 0.5 * (minX + maxX);
    torsionAxisZ = hingeZ;

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
     * 2. Изгиб:
     *
     * M_b*q_b'' + C_b*q_b' + K_b*q_b = Q_b
     */
    bendQDDot =
        (lastBendingForce
         - bendDamping * bendQDot
         - bendStiffness * bendQ)
        / bendMass;

    /*
     * 3. Кручение:
     *
     * I_t*q_t'' + C_t*q_t' + K_t*q_t = M_t
     */
    torsionQDDot =
        (lastTorsionMoment
         - torsionDamping * torsionQDot
         - torsionStiffness * torsionQ)
        / torsionInertia;

    /*
     * 4. Semi-implicit Euler.
     */
    bendQDot += bendQDDot * tau;
    bendQ    += bendQDot  * tau;

    torsionQDot += torsionQDDot * tau;
    torsionQ    += torsionQDot  * tau;

    /*
     * 5. Обновляем координаты узлов.
     */
    updateKinematics();

    /*
     * 6. Пересчитываем поля для записи в VTK.
     */
    updateAerodynamicFields();

    /*
     * 6. Пересчитываем напряжения в тетраэдрах
     */
    updateStressFields();

currentTime += tau;

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

void CalcMesh::configureWingTorsion(double initialTorsion,
                                    double initialTorsionVelocity,
                                    double torsionInertiaValue,
                                    double torsionStiffnessValue,
                                    double torsionDampingValue) {
    torsionQ = initialTorsion;
    torsionQDot = initialTorsionVelocity;
    torsionQDDot = 0.0;

    torsionInertia = std::max(torsionInertiaValue, 1.0e-9);
    torsionStiffness = std::max(torsionStiffnessValue, 0.0);
    torsionDamping = std::max(torsionDampingValue, 0.0);

    lastTorsionMoment = 0.0;

    updateKinematics();
    updateAerodynamicFields();
}

void CalcMesh::configureWingMaterial(double youngModulusValue,
                                     double poissonRatioValue,
                                     double yieldStressValue) {
    youngModulus = std::max(youngModulusValue, 1.0);
    poissonRatio = std::clamp(poissonRatioValue, 0.0, 0.49);

    shearModulus =
        youngModulus / (2.0 * (1.0 + poissonRatio));

    yieldStress = std::max(yieldStressValue, 1.0);

    maxEquivalentStress = 0.0;

    updateStressFields();
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

double CalcMesh::torsionShape(double relY) const {
    /*
     * relY = y0 - hingeY.
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
     * Линейная крутильная форма:
     *
     * psi(0) = 0  - у закрепления нет закрутки
     * psi(1) = 1  - torsionQ равен углу свободного конца
     *
     * Для первого шага линейная форма лучше кубической:
     * кручение обычно распределяется примерно линейно вдоль консоли.
     */
    return xi;
}

double CalcMesh::bendingShapeSecondDerivative(double relY) const {
    double xi = -relY / wingSpanLength;
    xi = std::clamp(xi, 0.0, 1.0);

    return (6.0 - 12.0 * xi) / (wingSpanLength * wingSpanLength);
}

double CalcMesh::torsionShapeDerivative(double relY) const {
    double xi = -relY / wingSpanLength;

    if (xi < 0.0 || xi > 1.0) {
        return 0.0;
    }

    return -1.0 / wingSpanLength;
}

void CalcMesh::updateStressFields() {
    elementEquivalentStress.assign(elements.size(), 0.0);
    maxEquivalentStress = 0.0;

    for (std::size_t elementIndex = 0; elementIndex < elements.size(); ++elementIndex) {
        const Element& element = elements[elementIndex];

        const CalcNode& n0 = nodes[element.nodesIds[0]];
        const CalcNode& n1 = nodes[element.nodesIds[1]];
        const CalcNode& n2 = nodes[element.nodesIds[2]];
        const CalcNode& n3 = nodes[element.nodesIds[3]];

        /*
         * Центр тетраэдра в текущей деформированной конфигурации.
         */
        const double centerX =
            0.25 * (n0.x + n1.x + n2.x + n3.x);

        const double centerY =
            0.25 * (n0.y + n1.y + n2.y + n3.y);

        const double centerZ =
            0.25 * (n0.z + n1.z + n2.z + n3.z);

        /*
         * Координата вдоль крыла относительно закрепления.
         */
        const double relY =
            centerY - hingeY;

        /*
         * Координаты относительно оси кручения.
         */
        const double relToTorsionAxisX =
            centerX - torsionAxisX;

        const double relToTorsionAxisZ =
            centerZ - torsionAxisZ;

        /*
         * Напряжение от изгиба:
         *
         * w(y) = bendQ * phi(y)
         * curvature = d2w/dy2 = bendQ * d2phi/dy2
         *
         * sigma = E * z * curvature
         */
        const double curvature =
            bendQ * bendingShapeSecondDerivative(relY);

        const double sigmaBending =
            youngModulus * relToTorsionAxisZ * curvature;

        /*
         * Напряжение от кручения:
         *
         * theta(y) = torsionQ * psi(y)
         * thetaDerivative = dtheta/dy
         *
         * tau = G * r * dtheta/dy
         */
        const double thetaDerivative =
            torsionQ * torsionShapeDerivative(relY);

        const double radiusFromTorsionAxis =
            std::sqrt(
                relToTorsionAxisX * relToTorsionAxisX +
                relToTorsionAxisZ * relToTorsionAxisZ
            );

        const double tauTorsion =
            shearModulus * radiusFromTorsionAxis * thetaDerivative;

        /*
         * Эквивалентное напряжение von Mises
         * для комбинации нормального и касательного напряжения.
         */
        const double equivalentStress =
            std::sqrt(
                sigmaBending * sigmaBending +
                3.0 * tauTorsion * tauTorsion
            );

        elementEquivalentStress[elementIndex] = equivalentStress;

        maxEquivalentStress =
            std::max(maxEquivalentStress, equivalentStress);
    }
}

void CalcMesh::updateKinematics() {
    com.centerX = 0.0;
    com.centerY = 0.0;
    com.centerZ = 0.0;

    for (auto& node : nodes) {
        const double phiB = bendingShape(node.relY);
        const double phiT = torsionShape(node.relY);

        /*
         * Исходные координаты узла.
         */
        const double x0 = hingeX + node.relX;
        const double y0 = hingeY + node.relY;
        const double z0 = hingeZ + node.relZ;

        /*
         * Изгиб по Z.
         */
        const double bendingDz = bendQ * phiB;
        const double bendingVz = bendQDot * phiB;

        /*
         * Локальный угол кручения данного сечения.
         */
        const double theta = torsionQ * phiT;
        const double thetaDot = torsionQDot * phiT;

        /*
         * Координаты относительно оси кручения.
         *
         * Кручение идет вокруг оси Y.
         */
        const double dx = x0 - torsionAxisX;
        const double dz = z0 - torsionAxisZ;

        /*
         * Малый поворот вокруг оси Y:
         *
         * x' = x + dz * theta
         * z' = z - dx * theta
         *
         * Плюс изгиб добавляется в z.
         */
        node.x = x0 + dz * theta;
        node.y = y0;
        node.z = z0 + bendingDz - dx * theta;

        /*
         * Скорости узлов:
         *
         * vx = dz * thetaDot
         * vz = bendingVz - dx * thetaDot
         */
        node.vx = dz * thetaDot;
        node.vy = 0.0;
        node.vz = bendingVz - dx * thetaDot;

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
    lastTorsionMoment = 0.0;

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

        // СДВИГ
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
        
        // КРУЧЕНИЕ
        const double phiB = bendingShape(node.relY);
        const double phiT = torsionShape(node.relY);

        /*
        * Обобщенная сила изгиба:
        *
        * z_i = z0_i + bendQ * phiB
        * dz_i / d(bendQ) = phiB
        */
        lastBendingForce += node.fz * phiB;

        /*
        * Обобщенный момент кручения.
        *
        * Для кручения вокруг оси Y:
        *
        * M_y = r_z * F_x - r_x * F_z
        *
        * где:
        *   r_x = x - torsionAxisX
        *   r_z = z - torsionAxisZ
        *
        * Так как локальный угол равен torsionQ * phiT,
        * обобщенный момент:
        *
        * Q_torsion = sum(M_y * phiT)
        */
        const double rx = node.x - torsionAxisX;
        const double rz = node.z - torsionAxisZ;

        const double momentY = rz * node.fx - rx * node.fz;

        lastTorsionMoment += momentY * phiT;
    }
}



void CalcMesh::snapshot(unsigned int snap_number) {
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    auto pressure = vtkSmartPointer<vtkDoubleArray>::New();
    pressure->SetName("pressure");
    
    auto stressArray = vtkSmartPointer<vtkDoubleArray>::New();
    stressArray->SetName("equivalent_stress");
    stressArray->SetNumberOfComponents(1);
    stressArray->SetNumberOfTuples(elements.size());

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

    for (std::size_t i = 0; i < elements.size(); ++i) {
        double stress = 0.0;

        if (i < elementEquivalentStress.size()) {
            stress = elementEquivalentStress[i];
        }

        stressArray->SetValue(i, stress);
    }

    unstructuredGrid->GetCellData()->AddArray(stressArray);

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
