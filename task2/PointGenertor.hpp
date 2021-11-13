#pragma once

#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <cfloat>

#include "Entities.hpp"


enum FillMode {
    Ascending,
    Descending,
    Random,
    Trigonometry
};

FillMode GetFillMode(const char* mode) {
    if (!strcmp("ascending", mode) || !strcmp("a", mode))
        return Ascending;
    if (!strcmp("descending", mode) || !strcmp("d", mode))
        return Descending;
    if (!strcmp("random", mode) || !strcmp("r", mode))
        return Random;
    if (!strcmp("trigonometry", mode) || !strcmp("t", mode))
        return Trigonometry;

    std::cout << "Warning: invalid fill mode (" << mode << "). Use trigonometry by default" << std::endl;
    return Trigonometry;
}

Point GetPoint(int i, int j, int n1, int n2, FillMode mode) {
    int total = n1 * n2;
    int index = i * n2 + j;
    float x, y;

    if (index >= total) {
        return (Point) {{FLT_MAX, FLT_MAX}, -1};
    }

    if (mode == Ascending) {
        x = index;
        y = index;
    }
    else if (mode == Descending) {
        x = total - index - 1;
        y = x;
    }
    else if (mode == Random) {
        x = rand() % total;
        y = rand() % total;
    }
    else {
        x = sin(index * 2 * M_PI / total);
        y = cos(index * 2 * M_PI / total);
    }

    return (Point) {{x, y}, index};
}
