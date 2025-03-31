#ifndef ELEMUNIV_H
#define ELEMUNIV_H

#include <vector>
#include "wezly.h"  
struct Surface {
    std::vector<std::vector<double>> N;

    Surface() {}

   
};
struct ElemUniv {
    std::vector<std::vector<double>> dN_dKsi;
    std::vector<std::vector<double>> dN_dEta;
    std::vector<Surface> surface;  
    std::vector<std::vector<double>> N;
    

    ElemUniv(int npc, const wezly& gauss);
    void print();
};

#endif // ELEMUNIV_H
