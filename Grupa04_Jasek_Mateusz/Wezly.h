#ifndef WEZLY_H
#define WEZLY_H

#include <vector>

struct wezly {
    std::vector<double> nodes;
    std::vector<double> weights;
    int DEGREE = 0;

    wezly(int degree);
};

#endif // WEZLY_H
