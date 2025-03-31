#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

struct node {
    double x;
    double y;
    int BC; 
    double Temp;
    node(double x, double y, int BC = 0) : x(x), y(y), BC(BC) {}

};

vector<node> NloadData(int nN);

#endif // NODE_H
