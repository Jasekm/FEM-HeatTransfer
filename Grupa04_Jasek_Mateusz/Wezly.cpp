#include "wezly.h"
#include <iostream>
#include <cmath>

using namespace std;  

wezly::wezly(int degree) {
    switch (degree) {
    case 0:
        nodes = { 0.0 };
        weights = { 2.0 };
        DEGREE = degree;
        break;
    case 1:
        nodes = { -1.0 / sqrt(3), 1.0 / sqrt(3) };
        weights = { 1.0, 1.0 };
        DEGREE = degree;
        break;
    case 2:
        nodes = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
        weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
        DEGREE = degree;
        break;
    case 3:
        nodes = { -0.861136, -0.339981, 0.339981, 0.861136 };
        weights = { 0.347855, 0.652145, 0.652145, 0.347855 };
        DEGREE = degree;
        break;
    case 4:
        nodes = { -0.906180, -0.538469, 0, 0.538469, 0.906180 };
        weights = { 0.236927, 0.478629, 0.568889, 0.478629, 0.236927 };
        DEGREE = degree;
        break;
    default:
        cout << "Stopieñ kwadratury nieobs³ugiwany" << endl;
    }
}
