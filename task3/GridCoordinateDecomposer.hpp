#pragma once
#include <vector>

#include <mpi.h>
#include "Entities.hpp"
#include "BatcherSorter.hpp"


class GridCoordinateDecomposer {
    int size;
    int rank;
    int n1, n2;
    MPI_Datatype pointType;

    void SetDomain(vector<Point> &points, int begin, int end, int domain);
    void Decompose(vector<Point> &points, int domain, int k, int begin, int n, int coord, DecomposeType type);
    int GetIntersections(vector<Point> &points, int begin, int leftCount, int rightCount);
public:
    GridCoordinateDecomposer(int size, int rank, int n1, int n2);

    void Decompose(vector<Point> &points, int k, DecomposeType type);
    vector<Point> JoinPoints(vector<Point>& points);
};

GridCoordinateDecomposer::GridCoordinateDecomposer(int size, int rank, int n1, int n2) {
    this->size = size;
    this->rank = rank;
    this->n1 = n1;
    this->n2 = n2;

    pointType = MakePointType();
}

void GridCoordinateDecomposer::SetDomain(vector<Point> &points, int begin, int end, int domain) {
    for (int i = 0; i < points.size(); i++) {
        int position = i + rank * points.size();

        if (position >= begin && position < end)
            points[i].domain = domain;
    }
}

// TODO
int GridCoordinateDecomposer::GetIntersections(vector<Point> &points, int begin, int leftCount, int rightCount) {
    return 0;
}

void GridCoordinateDecomposer::Decompose(vector<Point> &points, int k, DecomposeType type) {
    Decompose(points, 0, k, 0, n1 * n2, 0, type);
}

void GridCoordinateDecomposer::Decompose(vector<Point> &points, int domain, int k, int begin, int n, int coord, DecomposeType type) {
    if (k == 1) {
        SetDomain(points, begin, begin + n, domain);
        return;
    }

    int nextK = (k + 1) / 2;
    int nextN = n * ((double) nextK) / k;

    ParallelBetcherSorter sorter(rank, size, pointType);

    if (type == ChangeAxis) {
        sorter.Sort(points, PointComparator(coord, begin, begin + n));
        coord = !coord;
    }
    else if (type == Intersections) {
        vector<Point> tmpPoints(points);

        sorter.Sort(points, PointComparator(0, begin, begin + n));
        int countX = GetIntersections(points, begin, nextN, n - nextN);

        sorter.Sort(tmpPoints, PointComparator(1, begin, begin + n));
        int countY = GetIntersections(tmpPoints, begin, nextN, n - nextN);

        if (countY < countX) {
            points = tmpPoints;
        }
    }

    Decompose(points, domain, nextK, begin, nextN, coord, type);
    Decompose(points, domain + nextK, k - nextK, begin + nextN, n - nextN, coord, type);
}

vector<Point> GridCoordinateDecomposer::JoinPoints(vector<Point>& points) {
    int n = points.size();

    if (rank != 0) {
        MPI_Send(points.data(), n, pointType, 0, 0, MPI_COMM_WORLD);
        return points;
    }

    vector<Point> data(n * size);
    copy(points.begin(), points.end(), data.begin());

    for (int i = 1; i < size; i++) 
        MPI_Recv(data.data() + n * i, n, pointType, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return data;
}
