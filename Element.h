#ifndef ELEMENT_H
#define ELEMENT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "Jakobian.h"
class MacierzH;
using namespace std;
struct element {
    int ID[4];
    Jakobian jacobian;  
    MacierzH* macierzH;
    element(int a, int b, int c, int d) : ID{ a, b, c, d }, jacobian(), macierzH(nullptr) {}
};

vector<element> EloadData(int nE);

#endif // ELEMENT_H
