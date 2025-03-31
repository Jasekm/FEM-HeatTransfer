#pragma once
#ifndef JAKOBIAN_H
#define JAKOBIAN_H
#include <vector>
#include "eigen-master/Eigen/Dense"
#include "ElemUniv.h"
using namespace std;
struct grid;
struct Jakobian {
    vector<vector<double>> J;  
    vector<vector<double>> J1;
    vector<double> detJ;

    Jakobian();
    Jakobian(const ElemUniv& elem, const grid& grid, int pkt,int npc);
    void print();
};

#endif // JAKOBIAN_H
