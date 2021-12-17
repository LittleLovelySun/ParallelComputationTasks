#pragma once
#include <vector>

#include <mpi.h>
#include "Entities.hpp"
#include "BatcherSorter.hpp"


class GridCoordinateDecomposer {
    int size;
    int rank;
    MPI_Datatype pointType;

    void Decompose(vector<Point> &points, int domain, int k, int begin, int n, int coord);
    void FillPositions(vector<Point> &points);
public:
    GridCoordinateDecomposer(int size, int rank);
    void Decompose(vector<Point> &points, int k);
    vector<Point> JoinPoints(vector<Point>& points);
};

GridCoordinateDecomposer::GridCoordinateDecomposer(int size, int rank) {
    this->size = size;
    this->rank = rank;

    pointType = MakePointType();
}

void GridCoordinateDecomposer::FillPositions(vector<Point> &points) {
    for (size_t i = 0; i < points.size(); i++)
        points[i].position = i + rank * points.size();
}

void GridCoordinateDecomposer::Decompose(vector<Point> &points, int k) {
    FillPositions(points);
    Decompose(points, 0, k, 0, points.size() * size, 0);
}

void GridCoordinateDecomposer::Decompose(vector<Point> &points, int domain, int k, int begin, int n, int coord) {
    if (k == 1) {
        for (size_t i = 0; i < points.size(); i++) {
            int position = i + rank * points.size();
            if (position >= begin && position < begin + n) {
                points[i].domain = domain;
            }
        }

        return;
    }

    int k1 = (k + 1) / 2;
    int n1 = n * ((double) k1) / k;
    coord = !coord;

    ParallelBetcherSorter sorter(rank, size, pointType);
    FillPositions(points);

    PointComparator comparator(coord, begin, begin + n);
    sorter.Sort(points, comparator);

    Decompose(points, domain, k1, begin, n1, coord);
    Decompose(points, domain + k1, k - k1, begin + n1, n - n1, coord);
}

vector<Point> GridCoordinateDecomposer::JoinPoints(vector<Point>& points) {
    size_t n = points.size();

    if (rank != 0) {
        MPI_Send(points.data(), n, pointType, 0, 0, MPI_COMM_WORLD);
        return points;
    }

    size_t total = n * size;
    vector<Point> data(total);
    copy(points.begin(), points.end(), data.begin());

    for (int i = 1; i < size; i++) 
        MPI_Recv(data.data() + n * i, n, pointType, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return data;
}
