#pragma once

#include <algorithm>
#include <vector>
#include <mpi.h>

#include "Entities.hpp"


using namespace std;

struct PointComparator {
    int sortedCoord;
    int begin;
    int end;

    PointComparator(int sortedCoord, int begin, int end);
    bool operator()(const Point &a, const Point &b) const;
};

PointComparator::PointComparator(int sortedCoord, int begin, int end) {
    this->sortedCoord = sortedCoord;
    this->begin = begin;
    this->end = end;
}

bool PointComparator::operator()(const Point &p1, const Point &p2) const {
    if (p1.position >= begin && p1.position < end && p2.position >= begin && p2.position < end) {
        if (p1.coord[sortedCoord] != p2.coord[sortedCoord])
            return p1.coord[sortedCoord] < p2.coord[sortedCoord];
        
        return p1.coord[1 - sortedCoord] < p2.coord[1 - sortedCoord];
    }

    return p1.position < p2.position;
}


class ParallelBetcherSorter {
    int rank;
    int size;
    MPI_Datatype pointType;

    void MergePoints(vector<Point> &pointsI, int pi, int pj, const PointComparator &comparator);
    void Heapify(vector<Point> &points, int n, int i, const PointComparator &comparator);
    void HeapSort(vector<Point> &points, const PointComparator &comparator);
    void SetPositions(vector<Point> &points);
public:
    ParallelBetcherSorter(int rank, int size, MPI_Datatype pointType);
    void Sort(vector<Point> &points, const PointComparator &comparator);
    bool Check(vector<Point> &points, const PointComparator &comparator);
};

ParallelBetcherSorter::ParallelBetcherSorter(int rank, int size, MPI_Datatype pointType) {
    this->rank = rank;
    this->size = size;
    this->pointType = pointType;
}

void ParallelBetcherSorter::MergePoints(vector<Point> &pointsI, int pi, int pj, const PointComparator &comparator) {
    if (rank != pi && rank != pj)
        return;

    int n = pointsI.size();
    vector<Point> pointsJ(n);
    int other = rank == pi ? pj : pi;

    MPI_Request request;
    MPI_Status status;

    MPI_Isend(pointsI.data(), n, pointType, other, 0, MPI_COMM_WORLD, &request);
    MPI_Recv(pointsJ.data(), n, pointType, other, 0, MPI_COMM_WORLD, &status);
    MPI_Wait(&request, &status);

    vector<Point> points(n);

    if (rank == pi) {
        for (int i = 0, j = 0, k = 0; k < n; k++) {
            points[k] = comparator(pointsI[i], pointsJ[j]) ? pointsI[i++] : pointsJ[j++];
        }
    } else {
        for (int i = n - 1, j = n - 1, k = n - 1; k >= 0; k--) {
            points[k] = comparator(pointsJ[j], pointsI[i]) ? pointsI[i--] : pointsJ[j--];
        }
    }

    pointsI = points;
}

bool ParallelBetcherSorter::Check(vector<Point> &points, const PointComparator &comparator) {
    int n = points.size();

    if (rank != 0) {
        MPI_Send(points.data(), n, pointType, 0, 0, MPI_COMM_WORLD);
        return true;
    }

    int total = n * size;
    vector<Point> data(total);
    copy(points.begin(), points.end(), data.begin());

    for (int i = 1; i < size; i++) 
        MPI_Recv(data.data() + n * i, n, pointType, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < total - 1; i++)
        if (comparator(data[i + 1], data[i]))
            return false;

    return true;
}

void ParallelBetcherSorter::Heapify(vector<Point> &points, int n, int i, const PointComparator &comparator) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < n && comparator(points[largest], points[left]))
        largest = left;

    if (right < n && comparator(points[largest], points[right]))
        largest = right;

    if (largest != i) {
        swap(points[i], points[largest]);
        Heapify(points, n, largest, comparator);
    }
}

void ParallelBetcherSorter::HeapSort(vector<Point> &points, const PointComparator &comparator) {
    for (int i = points.size() / 2 - 1; i >= 0; i--)
        Heapify(points, points.size(), i, comparator);

    for (int i = points.size() - 1; i > 0; i--) {
        swap(points[0], points[i]);
        Heapify(points, i, 0, comparator);
    }
}

void ParallelBetcherSorter::SetPositions(vector<Point> &points) {
    for (int i = 0; i < points.size(); i++)
        points[i].position = i + rank * points.size();
}

void ParallelBetcherSorter::Sort(vector<Point> &points, const PointComparator &comparator) {
    int n = points.size();
    int startRank = comparator.begin / n;
    int endRank = (comparator.end - 1) / n;

    SetPositions(points);

    if (startRank <= rank && rank <= endRank)
        HeapSort(points, comparator);

    for (int p = 1; p < size; p *= 2)
        for (int k = p; k >= 1; k /= 2)
            for (int j = k % p; j <= size - 1 - k; j += 2*k)
                for (int i = 0; i < min(k, size - k - j); i++)
                    if ((i + j) / (p * 2) == (i + j + k) / (p * 2))
                        if ((startRank <= i + j && i + j <= endRank) || (startRank <= i + j + k && i + j + k <= endRank))
                            MergePoints(points, i + j, i + j + k, comparator);
}
