#pragma once

#include <algorithm>
#include <vector>
#include <mpi.h>

#include "Entities.hpp"


using namespace std;

class PointComparator {
    int sortedCoord;
    int begin;
    int end;
public:
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

    void MakeMergePoints(const vector<Point> &pointsI, const vector<Point> &pointsJ, vector<Point> &points, const PointComparator &comparator);
    void MergePoints(vector<Point> &pointsI, int pi, int pj, const PointComparator &comparator);
    void Merge(int leftPos, int rightPos, int leftSize, int rightSize, int offset, vector<Point> &points, const PointComparator &comparator);

    void Sort(int position, int length, vector<Point> &points, const PointComparator &comparator);
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

void ParallelBetcherSorter::MakeMergePoints(const vector<Point> &pointsI, const vector<Point> &pointsJ, vector<Point> &points, const PointComparator &comparator) {
    size_t i = 0;
    size_t j = 0;
    size_t k = 0;

    while (i < pointsI.size() && j < pointsJ.size())
        points[k++] = comparator(pointsI[i], pointsJ[j]) ? pointsI[i++] : pointsJ[j++];

    while (i < pointsI.size())
        points[k++] = pointsI[i++];

    while (j < pointsJ.size())
        points[k++] = pointsJ[j++];
}

void ParallelBetcherSorter::MergePoints(vector<Point> &pointsI, int pi, int pj, const PointComparator &comparator) {
    size_t n = pointsI.size();
    
    if (rank == pj) {
        MPI_Send(pointsI.data(), n, pointType, pi, 0, MPI_COMM_WORLD);
        MPI_Recv(pointsI.data(), n, pointType, pi, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        return;
    }

    if (pi != rank)
        return;

    vector<Point> pointsJ(n);
    MPI_Recv(pointsJ.data(), n, pointType, pj, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    vector<Point> points(n * 2);

    MakeMergePoints(pointsI, pointsJ, points, comparator);
    MPI_Send(points.data() + n, n, pointType, pj, 0, MPI_COMM_WORLD);

    for (int i = 0; i < n; i++)
        pointsI[i] = points[i];
}

void ParallelBetcherSorter::Merge(int leftPos, int rightPos, int leftSize, int rightSize, int offset, vector<Point> &points, const PointComparator &comparator) {
	if (leftSize == 0 || rightSize == 0)
		return;

	if (leftSize == 1 && rightSize == 1) {
		MergePoints(points, leftPos, rightPos, comparator);
		return;
	}

	Merge(leftPos, rightPos, (leftSize + 1) / 2, (rightSize + 1) / 2, offset * 2, points, comparator);
	Merge(leftPos + offset, rightPos + offset, leftSize / 2, rightSize / 2, offset * 2, points, comparator);

	for (int i = 1; i < leftSize - 1; i += 2)
		MergePoints(points, leftPos + i * offset, leftPos + (i + 1) * offset, comparator);

	if (leftSize % 2 == 0)
		MergePoints(points, leftPos + (leftSize - 1) * offset, rightPos, comparator);

	for (int i = 1 - leftSize % 2; i < rightSize - 1; i += 2)
		MergePoints(points, rightPos + i * offset, rightPos + (i + 1) * offset, comparator);
}

void ParallelBetcherSorter::Sort(int position, int length, vector<Point> &points, const PointComparator &comparator) {
	if (length < 2)
		return;

	int middle = length / 2;
	Sort(position, middle, points, comparator);
	Sort(position + middle, length - middle, points, comparator);
	Merge(position, position + middle, middle, length - middle, 1, points, comparator);
}

bool ParallelBetcherSorter::Check(vector<Point> &points, const PointComparator &comparator) {
    size_t n = points.size();

    if (rank != 0) {
        MPI_Send(points.data(), n, pointType, 0, 0, MPI_COMM_WORLD);
        return true;
    }

    size_t total = n * size;
    vector<Point> data(total);
    copy(points.begin(), points.end(), data.begin());

    for (int i = 1; i < size; i++) 
        MPI_Recv(data.data() + n * i, n, pointType, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < total - 1; i++)
        if (comparator(data[i + 1], data[i]))
            return false;

    return true;
}

void ParallelBetcherSorter::Sort(vector<Point> &points, const PointComparator &comparator) {
    sort(points.begin(), points.end(), comparator);
    Sort(0, size, points, comparator);
}
