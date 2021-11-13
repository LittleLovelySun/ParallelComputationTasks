#pragma once

#include <iostream>
#include <vector>
#include <cstring>
#include "Entities.hpp"


enum FillMode {
    Ascending,
    Descending,
    Random
};

FillMode GetFillMode(const char* mode) {
    if (!strcmp("ascending", mode) || !strcmp("a", mode))
        return Ascending;
    if (!strcmp("descending", mode) || !strcmp("d", mode))
        return Descending;
    if (!strcmp("random", mode) || !strcmp("r", mode))
        return Random;

    std::cout << "Warning: invalid fill mode (" << mode << "). Use random by default" << std::endl;
    return Random;
}

Point GetPoint(int i, int j, int n1, int n2, FillMode mode) {
    int index = i * n2 + j;
    float x, y;
    if (mode == Ascending) {
        x = index;
        y = index;
    }
    else if (mode == Descending) {
        x = n1 * n2 - index - 1;
        y = x;
    }
    else {
        x = rand() % (n1 * n2);
        y = rand() % (n1 * n2);
    }

    return (Point) {{x, y}, index};
}
