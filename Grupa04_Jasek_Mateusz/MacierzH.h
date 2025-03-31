// MacierzH.h
#ifndef MACIERZH_H
#define MACIERZH_H

#include <vector>
#include <iostream>
#include "ElemUniv.h"
#include "Jakobian.h"
#include "GlobalData.h"
#include "Grid.h"
struct grid;
using namespace std;
class MacierzH {
public:
    vector<vector<double>> macierz;
    vector<vector<double>> dN_dx;
    vector<vector<double>> dN_dy;
    vector<vector<double>> macierzHbc;
    vector<vector<double>> macierzC;
    vector<double> P;

    MacierzH();
    MacierzH(const ElemUniv& elem, const Jakobian& jacobian, int npc, const wezly& gauss, GlobalData& data,const grid& grid,int elemNumber);

    void print_dN_dx() const;
    void print_dN_dy() const;
    void print_macierzH() const;
    void print_macierzC() const;

private:
    vector<vector<double>> transpose(const vector<double>& vec);
    vector<vector<double>> multiply(const vector<double>& vec, const vector<vector<double>>& transposed_vec);
    vector<vector<double>> add(const vector<vector<double>>& A, const vector<vector<double>>& B);
    vector<vector<double>> multiplyByScalar(const vector<vector<double>>& matrix, double scalar);
    vector<double> multiplyVectorByScalar(const std::vector<double>& vec, double number);
    vector<double> addVectors(const std::vector<double>& vec1, const std::vector<double>& vec2);
    void print(const vector<vector<double>>& matrix);
    void printVector(const std::vector<double>& vec);
};

#endif // MACIERZH_H
