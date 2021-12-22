#pragma once
#include <vector>

#include <mpi.h>
#include "Entities.hpp"
#include "BatcherSorter.hpp"


class GridCoordinateDecomposer {
    int size;
    int rank;
    int n1, n2;
    int totalEdges;
    MPI_Datatype pointType;

    void SetDomain(vector<Point> &points, int begin, int end, int domain);
    void Decompose(vector<Point> &points, int domain, int k, int begin, int n, int coord, DecomposeType type);

    vector<Point> JoinPoints(vector<Point>& points);
    void FindTotalEdges(vector<Point> &points);

    int GetSingleIntersections(vector<Point> &points, int begin, int leftCount, int rightCount);
    int GetIntersections(vector<Point> &points, int begin, int leftCount, int rightCount);
public:
    GridCoordinateDecomposer(int size, int rank, int n1, int n2);

    void Decompose(vector<Point> &points, int k, DecomposeType type);
    int GetTotalEdges(vector<Point>& points, DecomposeType type);
    void Save(vector<Point>& points, const char *path);
};

GridCoordinateDecomposer::GridCoordinateDecomposer(int size, int rank, int n1, int n2) {
    this->size = size;
    this->rank = rank;
    this->n1 = n1;
    this->n2 = n2;

    pointType = MakePointType();
    totalEdges = 0;
}

void GridCoordinateDecomposer::SetDomain(vector<Point> &points, int begin, int end, int domain) {
    for (int i = 0; i < points.size(); i++) {
        int position = i + rank * points.size();

        if (position >= begin && position < end)
            points[i].domain = domain;
    }
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

void GridCoordinateDecomposer::FindTotalEdges(vector<Point> &points) {
    vector<Point> totalPoints = JoinPoints(points);

    totalEdges = 0;

    for (int i = 0; i < totalPoints.size(); i++)
        for (int j = i + 1; j < totalPoints.size(); j++)
            if (Distance(totalPoints[i], totalPoints[j], n2) == 1 && totalPoints[i].domain != totalPoints[j].domain)
                totalEdges++;

    MPI_Bcast(&totalEdges, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

int GridCoordinateDecomposer::GetSingleIntersections(vector<Point> &points, int begin, int leftCount, int rightCount) {
    int count = 0;

    for (int i = 0; i < leftCount; i++)
        for (int j = 0; j < rightCount; j++)
            if (Distance(points[begin + i], points[begin + leftCount + j], n2) == 1)
                count++;

    return count;
}

int GridCoordinateDecomposer::GetIntersections(vector<Point> &points, int begin, int leftCount, int rightCount) {
    int n = points.size();

    int beginRank = begin / n;
    int endRank = (begin + leftCount + rightCount - 1) / n;

    // если весь интервал внезапно на одном процессе
    if (beginRank == endRank) {
        int count = 0;

        if (rank == beginRank)
            count = GetSingleIntersections(points, begin % n, leftCount, rightCount);

        MPI_Bcast(&count, 1, MPI_INT, beginRank, MPI_COMM_WORLD);
        return count;
    }
    
    int endLeftRank = (begin + leftCount - 1) / n;
    int beginRightRank = (begin + leftCount) / n;
    int count = 0;

    if (beginRightRank <= rank && rank <= endRank) {
        MPI_Request request;

        for (int i = beginRank; i <= endLeftRank; i++)
            MPI_Isend(points.data(), n, pointType, i, 0, MPI_COMM_WORLD, &request);
    }
    
    if (beginRank <= rank && rank <= endLeftRank) {
        int total = n * (endRank - beginRightRank + 1);

        vector<Point> rightPoints(total);
        vector<MPI_Request> requests(endRank - beginRightRank + 1);
        vector<MPI_Status> statuses(endRank - beginRightRank + 1);

        for (int i = beginRightRank; i <= endRank; i++)
            MPI_Irecv(rightPoints.data() + (i - beginRightRank) * n, n, pointType, i, 0, MPI_COMM_WORLD, &requests[i - beginRightRank]);
        MPI_Waitall(requests.size(), requests.data(), statuses.data());

        int beginI = rank == beginRank ? begin % n : 0;
        int endI = rank == endLeftRank ? (begin + leftCount - 1) % n : n - 1;

        for (int i = beginI; i <= endI; i++) {
            for (int rankJ = beginRightRank; rankJ <= endRank; rankJ++) {
                int beginJ = rankJ == beginRightRank ? (begin + leftCount) % n : 0;
                int endJ = rankJ == endRank ? (begin + leftCount + rightCount - 1) % n : n - 1;

                for (int j = beginJ; j <= endJ; j++)
                    if (Distance(points[i], rightPoints[j + n * (rankJ - beginRightRank)], n2) == 1)
                        count++;
            }
        }
    }

    int totalCount = 0;
    MPI_Reduce(&count, &totalCount, 1, MPI_INT, MPI_SUM, beginRank, MPI_COMM_WORLD);
    MPI_Bcast(&totalCount, 1, MPI_INT, beginRank, MPI_COMM_WORLD);

    return totalCount;
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
            totalEdges += countY;
        }
        else {
            totalEdges += countX;
        }
    }

    Decompose(points, domain, nextK, begin, nextN, coord, type);
    Decompose(points, domain + nextK, k - nextK, begin + nextN, n - nextN, coord, type);
}

int GridCoordinateDecomposer::GetTotalEdges(vector<Point>& points, DecomposeType type) {
    if (type == ChangeAxis)
        FindTotalEdges(points);

    return totalEdges;
}

void GridCoordinateDecomposer::Save(vector<Point>& points, const char *path) {
    for (int p = 0; p < size; p++) {
        if (rank == p) {
            ofstream fout(path, rank == 0 ? ios::out : ios::app);
            PrintPoints(points, n2, fout);
            fout.close();
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}
