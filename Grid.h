#ifndef GRID_H
#define GRID_H

#include <vector>
#include <iostream>
#include "Element.h"
#include "Node.h"

using namespace std;

struct grid {
    int nN;
    int nE;
    vector<element> elements;
    vector<node> nodes;

    grid(int nNodes, int nElements, const vector<element>& elems, const vector<node>& nods)
        : nN(nNodes), nE(nElements), elements(elems), nodes(nods) {}
};

void displayGrid(grid newGrid);

#endif // GRID_H
