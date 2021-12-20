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
