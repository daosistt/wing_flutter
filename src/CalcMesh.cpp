#include "CalcMesh.h"
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <cmath>
#include <string>

CalcMesh::CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints) : currentTime(0.0)  {
    nodes.resize(nodesCoords.size() / 3);
    for(unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
        double pointX = nodesCoords[i*3];
        double pointY = nodesCoords[i*3 + 1];
        double pointZ = nodesCoords[i*3 + 2];

        // TODO
        double smth = pow(pointX, 2) + pow(pointY, 2) + pow(pointZ, 2);
        nodes[i] = CalcNode(pointX, pointY, pointZ, smth, 0.0, 0.0, 0.0);
    }
    
    updateCenterOfMass(); 

    for (auto& node : nodes) {
        node.relX = node.x - com.centerX;
        node.relY = node.y - com.centerY;
        node.relZ = node.z - com.centerZ;
    }

    elements.resize(tetrsPoints.size() / 4);
    for(unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
        elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
        elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
        elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
        elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
    }
}




void CalcMesh::doTimeStep(double tau) {
    for(auto& node : nodes) {
        node.vx = vxFunc(node.x, node.y, node.z, currentTime);
        node.vy = vyFunc(node.x, node.y, node.z, currentTime);
        node.vz = vzFunc(node.x, node.y, node.z, currentTime);
    }

    for(auto& node : nodes) {
        node.move(tau);
    }

    updateCenterOfMass(); 
    com.yaw   += com.yaw_rate_func(currentTime) * tau;
    com.pitch += com.pitch_rate_func(currentTime) * tau;
    com.roll  += com.roll_rate_func(currentTime) * tau;

    double cosYaw = cos(com.yaw);
    double sinYaw = sin(com.yaw);
    double cosPitch = cos(com.pitch);
    double sinPitch = sin(com.pitch);
    double cosRoll = cos(com.roll);
    double sinRoll = sin(com.roll);
    for (auto& node : nodes) {
        double x = node.relX;
        double y = node.relY;
        double z = node.relZ;

        double x1 = x * cosYaw - y * sinYaw;
        double y1 = x * sinYaw + y * cosYaw;
        double z1 = z;

        double x2 = x1 * cosPitch + z1 * sinPitch;
        double y2 = y1;
        double z2 = -x1 * sinPitch + z1 * cosPitch;

        double x3 = x2;
        double y3 = y2 * cosRoll - z2 * sinRoll;
        double z3 = y2 * sinRoll + z2 * cosRoll;

        node.x = com.centerX + x3;
        node.y = com.centerY + y3;
        node.z = com.centerZ + z3;
    }


    currentTime += tau;
}



void CalcMesh::snapshot(unsigned int snap_number) {
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    auto smth = vtkSmartPointer<vtkDoubleArray>::New();
    smth->SetName("smth");

    auto vel = vtkSmartPointer<vtkDoubleArray>::New();
    vel->SetName("velocity");
    vel->SetNumberOfComponents(3);

    for(unsigned int i = 0; i < nodes.size(); i++) {
        dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);
        double _vel[3] = {nodes[i].vx, nodes[i].vy, nodes[i].vz};
        vel->InsertNextTuple(_vel);
        smth->InsertNextValue(nodes[i].smth);
    }

    unstructuredGrid->SetPoints(dumpPoints);
    unstructuredGrid->GetPointData()->AddArray(vel);
    unstructuredGrid->GetPointData()->AddArray(smth);

    for(unsigned int i = 0; i < elements.size(); i++) {
        auto tetra = vtkSmartPointer<vtkTetra>::New();
        tetra->GetPointIds()->SetId(0, elements[i].nodesIds[0]);
        tetra->GetPointIds()->SetId(1, elements[i].nodesIds[1]);
        tetra->GetPointIds()->SetId(2, elements[i].nodesIds[2]);
        tetra->GetPointIds()->SetId(3, elements[i].nodesIds[3]);
        unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
    }

    std::string fileName = "jet-step-" + std::to_string(snap_number) + ".vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}



void CalcMesh::setVelocities(std::function<double(double,double,double,double)> fx,
                             std::function<double(double,double,double,double)> fy,
                             std::function<double(double,double,double,double)> fz) {
    for(auto& node : nodes) {
        vxFunc = fx;
        vyFunc = fy;
        vzFunc = fz;
    }
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

void CalcMesh::setAngularVelocity(
    std::function<double(double)> yaw_rate,
    std::function<double(double)> pitch_rate,
    std::function<double(double)> roll_rate) {

    com.yaw_rate_func = yaw_rate;
    com.pitch_rate_func = pitch_rate;
    com.roll_rate_func = roll_rate;
}
