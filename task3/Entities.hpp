#pragma once
#include <iostream>
#include <vector>
#include <cstdlib>


struct Point {
    float coord[2];
    int index;
    int position;
    int domain;
};

enum GridType {
    Grid,
    Spiral,
    Circle,
    Flower
};

enum DecomposeType {
    ChangeAxis,
    Intersections
};


void PrintPoints(const std::vector<Point> &points, int n2, std::ostream& out) {
    for (size_t index = 0; index < points.size(); index++) {
        if (points[index].index == -1)
            continue;

        int i = points[index].index / n2;
        int j = points[index].index % n2;
        float x = points[index].coord[0];
        float y = points[index].coord[1];

        out << i << " " << j << " " << x << " " << y << " " << points[index].domain << std::endl;
    }
}

int Distance(Point p1, Point p2, int n2) {
    int i1 = p1.index / n2;
    int j1 = p1.index % n2;

    int i2 = p2.index / n2;
    int j2 = p2.index % n2;

    return abs(i1 - i2) + abs(j1 - j2);
}

MPI_Datatype MakePointType() {
    Point point;
    int lengths[] = {2, 1, 1, 1};

    MPI_Aint displacements[4];
    MPI_Aint address;

    MPI_Get_address(&point, &address);
    MPI_Get_address(&point.coord[0], &displacements[0]);
    MPI_Get_address(&point.index, &displacements[1]);
    MPI_Get_address(&point.position, &displacements[2]);
    MPI_Get_address(&point.domain, &displacements[3]);

    for (int i = 0; i < 4; i++) 
        displacements[i] = displacements[i] - address;

    MPI_Datatype types[4] = {MPI_FLOAT, MPI_INT, MPI_INT, MPI_INT};

    MPI_Datatype pointType;
    MPI_Type_create_struct(4, lengths, displacements, types, &pointType);
    MPI_Type_commit(&pointType);

    return pointType;
}

GridType GetGridType(const char* type) {
    if (!strcmp("grid", type) || !strcmp("g", type))
        return Grid;

    if (!strcmp("spiral", type) || !strcmp("s", type))
        return Spiral;

    if (!strcmp("circle", type) || !strcmp("c", type))
        return Circle;

    if (!strcmp("flower", type) || !strcmp("f", type))
        return Flower;

    std::cout << "Warning: invalid fill type (" << type << "). Used grid by default" << std::endl;
    return Grid;
}

std::ostream& operator<<(std::ostream& os, GridType type) {
    if (type == Grid)
        os << "grid";
    else if (type == Spiral)
        os << "spiral";
    else if (type == Circle)
        os << "circle";
    else if (type == Flower)
        os << "flower";

    return os;
}

DecomposeType GetDecomposeType(const char* type) {
    if (!strcmp("change-axis", type) || !strcmp("c", type))
        return ChangeAxis;

    if (!strcmp("intersections", type) || !strcmp("i", type))
        return Intersections;

    std::cout << "Warning: invalid decompose type (" << type << "). Used intersections by default" << std::endl;
    return Intersections;
}

std::ostream& operator<<(std::ostream& os, DecomposeType type) {
    if (type == ChangeAxis)
        os << "change axis";
    else if (type == Intersections)
        os << "intersections";

    return os;
}