#pragma once

#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <cfloat>

#include "Entities.hpp"

using namespace std;

enum GrigType {
    Grid,
    Spiral,
    Circle,
    Flower
};


GrigType GetGridType(const char* type) {
    if (!strcmp("grid", type) || !strcmp("g", type))
        return Grid;
    if (!strcmp("spiral", type) || !strcmp("s", type))
        return Spiral;
    if (!strcmp("circle", type) || !strcmp("c", type))
        return Circle;
    if (!strcmp("flower", type) || !strcmp("f", type))
        return Flower;

    cout << "Warning: invalid fill type (" << type << "). Used grid by default" << endl;
    return Grid;
}


class GridGenerator {
    int size;
    int rank;
    GrigType type;

    Point GetPoint(int i, int j, int n1, int n2);
public:
    GridGenerator(int size, int rank, GrigType type);
    vector<Point> Generate(int n1, int n2);
};

GridGenerator::GridGenerator(int size, int rank, GrigType type) {
    this->size = size;
    this->rank = rank;
    this->type = type;
}

Point GridGenerator::GetPoint(int i, int j, int n1, int n2) {
    int index = i * n2 + j;
    float x, y;
    
    if (type == Grid) {
        x = i;
        y = j;
    }
    else if (type == Circle) {
        double t = 2 * M_PI * i / n1;
        double r = (j + 1.0) / n2;

        x = r * cos(t);
        y = r * sin(t);
    }
    else if (type == Spiral) {
        double t = -(i + 1.0) / n1 * M_PI * 6;
        double r = t - j * 2.0 / n2 - 1;

        x = r * cos(t);
        y = r * sin(t);
    }
    else {
        double t = i * 2.0 / n1 * M_PI;
        double r = cos(4 * t) * (j + 1.0) / n2;

        x = r * cos(t);
        y = r * sin(t);
    }

    return (Point) {{x, y}, index, -1, -1};
}

vector<Point> GridGenerator::Generate(int n1, int n2) {
    int total = n1 * n2;
    int count = (total + size - 1) / size;
    vector<Point> points(count);

    for (int index = 0; index < count; index++) {
        int position = index + rank * count;
        int i = position / n2;
        int j = position % n2;
        
        if (position < total)
            points[index] = GetPoint(i, j, n1, n2);
        else
            points[index] = (Point) {{FLT_MAX, FLT_MAX}, -1, -1, -1};
    }

    return points;
}