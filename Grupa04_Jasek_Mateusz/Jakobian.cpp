#include "Jakobian.h"
#include "grid.h"  
#include "ElemUniv.h"
using namespace std;  

Jakobian::Jakobian() : detJ() {}

Jakobian::Jakobian(const ElemUniv& elem, const grid& grid, int elemet_number, int npc) {
    J.resize(npc, vector<double>(4));  
    J1.resize(npc, vector<double>(4)); 
    detJ.resize(npc);

    for (int i = 0; i < npc; i++) {

        double dx_dKsi =
            elem.dN_dKsi[i][0] * grid.nodes[grid.elements[elemet_number].ID[0] - 1].x +
            elem.dN_dKsi[i][1] * grid.nodes[grid.elements[elemet_number].ID[1] - 1].x +
            elem.dN_dKsi[i][2] * grid.nodes[grid.elements[elemet_number].ID[2] - 1].x +
            elem.dN_dKsi[i][3] * grid.nodes[grid.elements[elemet_number].ID[3] - 1].x;

        double dx_dEta =
            elem.dN_dEta[i][0] * grid.nodes[grid.elements[elemet_number].ID[0] - 1].x +
            elem.dN_dEta[i][1] * grid.nodes[grid.elements[elemet_number].ID[1] - 1].x +
            elem.dN_dEta[i][2] * grid.nodes[grid.elements[elemet_number].ID[2] - 1].x +
            elem.dN_dEta[i][3] * grid.nodes[grid.elements[elemet_number].ID[3] - 1].x;

        double dy_dKsi =
            elem.dN_dKsi[i][0] * grid.nodes[grid.elements[elemet_number].ID[0] - 1].y +
            elem.dN_dKsi[i][1] * grid.nodes[grid.elements[elemet_number].ID[1] - 1].y +
            elem.dN_dKsi[i][2] * grid.nodes[grid.elements[elemet_number].ID[2] - 1].y +
            elem.dN_dKsi[i][3] * grid.nodes[grid.elements[elemet_number].ID[3] - 1].y;

        double dy_dEta =
            elem.dN_dEta[i][0] * grid.nodes[grid.elements[elemet_number].ID[0] - 1].y +
            elem.dN_dEta[i][1] * grid.nodes[grid.elements[elemet_number].ID[1] - 1].y +
            elem.dN_dEta[i][2] * grid.nodes[grid.elements[elemet_number].ID[2] - 1].y +
            elem.dN_dEta[i][3] * grid.nodes[grid.elements[elemet_number].ID[3] - 1].y;

        J[i][0] = dx_dKsi;
        J[i][1] = dy_dKsi;
        J[i][2] = dx_dEta;
        J[i][3] = dy_dEta;

        Eigen::Matrix2d J_matrix;
        J_matrix << J[i][0], J[i][1],
            J[i][2], J[i][3];

        detJ[i] = J_matrix.determinant();
        Eigen::Matrix2d J1_matrix = J_matrix.inverse();


        J1[i][0] = J1_matrix(0, 0);
        J1[i][1] = J1_matrix(0, 1);
        J1[i][2] = J1_matrix(1, 0);
        J1[i][3] = J1_matrix(1, 1);
    }
}

void Jakobian::print() {
    cout << "\nJakobian:" << endl;
    for (const auto& row : J) {
        cout << "[" << row[0] << ", " << row[1] << ", " << row[2] << ", " << row[3] << "]" << endl;
    }

    cout << "\nJakobian odwrotny:" << endl;
    for (const auto& row : J1) {
        cout << "[" << row[0] << ", " << row[1] << ", " << row[2] << ", " << row[3] << "]" << endl;
    }

    cout << "\nWyznaczniki macierzy:" << endl;
    for (const auto& determinant : detJ) {
        cout << determinant << " ";
    }
    cout << endl;
    cout << endl;
}