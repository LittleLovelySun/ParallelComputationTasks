#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <cstring>

#include <mpi.h>

#include "Entities.hpp"
#include "PointGenertor.hpp"

using namespace std;


class PointComparator {
    int sortedCoord;
public:
    PointComparator(int sortedCoord);
    bool operator()(const Point &a, const Point &b);
};

PointComparator::PointComparator(int sortedCoord) {
    this->sortedCoord = sortedCoord;
}

bool PointComparator::operator()(const Point &a, const Point &b) {
    return a.coord[sortedCoord] < b.coord[sortedCoord];
}


class ParallelBetcherSorter {
    int rank;
    int size;
    int sortedCoord;
    MPI_Datatype pointType;

    PointComparator comparator;

    void Heapify(vector<Point> &points, size_t n, size_t i);
    void HeapSort(vector<Point> &points);
    void MergeSort(vector<Point> &points, int left, int right);

    void MakePointType();
    void MakeMergePoints(const vector<Point> &pointsI, const vector<Point> &pointsJ, vector<Point> &points);
    void MergePoints(vector<Point> &pointsI, int pi, int pj);
    void Merge(int leftPos, int rightPos, int leftSize, int rightSize, int offset, vector<Point> &points);

    void Sort(int position, int length, vector<Point> &points);
public:
    ParallelBetcherSorter(int rank, int size, int sortedCoord);
    void Sort(vector<Point> &points, SortAlgorithm algorithm);
    bool Check(vector<Point> &points);
};

ParallelBetcherSorter::ParallelBetcherSorter(int rank, int size, int sortedCoord): comparator(sortedCoord) {
    this->rank = rank;
    this->size = size;
    this->sortedCoord = sortedCoord;

    MakePointType();
}

void ParallelBetcherSorter::MakePointType() {
    Point point;
    int lengths[] = {2, 1};

    MPI_Aint displacements[2];
    MPI_Aint address;

    MPI_Get_address(&point, &address);
    MPI_Get_address(&point.coord[0], &displacements[0]);
    MPI_Get_address(&point.index, &displacements[1]);

    displacements[0] = displacements[0] - address;
    displacements[1] = displacements[1] - address;

    MPI_Datatype types[2] = {MPI_FLOAT, MPI_INT};

    MPI_Type_create_struct(2, lengths, displacements, types, &pointType);
    MPI_Type_commit(&pointType);
}

void ParallelBetcherSorter::Heapify(vector<Point> &points, size_t n, size_t i) {
    size_t largest = i;

    while (true) {
        size_t left = 2 * i + 1;
        size_t right = 2 * i + 2;

        if (left < n && comparator(points[largest], points[left]))
            largest = left;

        if (right < n && comparator(points[largest], points[right]))
            largest = right;

        if (largest == i)
            return;

        swap(points[i], points[largest]);
        i = largest;
    }
}

void ParallelBetcherSorter::HeapSort(vector<Point> &points) {
    for (size_t i = points.size() / 2 - 1; i + 1 > 0; i--) {
        Heapify(points, points.size(), i);
    }

    for (size_t i = points.size() - 1; i + 1 > 0; i--) {
        swap(points[0], points[i]);
        Heapify(points, i, 0);
    }
}

void ParallelBetcherSorter::MergeSort(vector<Point> &points, int left, int right) {
    if (left + 1 > right)
        return;

    int mid = left + (right - left) / 2;
    MergeSort(points, left, mid);
    MergeSort(points, mid + 1, right);
    
    vector<Point> tmp(right - left + 1);

    size_t i = left;
    size_t j = mid + 1;
    size_t k = 0;

    while (i <= mid && j <= right)
        tmp[k++] = comparator(points[i], points[j]) ? points[i++] : points[j++];

    while (i <= mid)
        tmp[k++] = points[i++];

    while (j <= right)
        tmp[k++] = points[j++];

    for (size_t i = left; i <= right; i++)
        points[i] = tmp[i - left];
}

void ParallelBetcherSorter::MakeMergePoints(const vector<Point> &pointsI, const vector<Point> &pointsJ, vector<Point> &points) {
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

void ParallelBetcherSorter::MergePoints(vector<Point> &pointsI, int pi, int pj) {
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

    MakeMergePoints(pointsI, pointsJ, points);
    MPI_Send(points.data() + n, n, pointType, pj, 0, MPI_COMM_WORLD);

    for (int i = 0; i < n; i++)
        pointsI[i] = points[i];
}

void ParallelBetcherSorter::Merge(int leftPos, int rightPos, int leftSize, int rightSize, int offset, vector<Point> &points) {
	if (leftSize == 0 || rightSize == 0)
		return;

	if (leftSize == 1 && rightSize == 1) {
		MergePoints(points, leftPos, rightPos);
		return;
	}

	Merge(leftPos, rightPos, (leftSize + 1) / 2, (rightSize + 1) / 2, offset * 2, points);
	Merge(leftPos + offset, rightPos + offset, leftSize / 2, rightSize / 2, offset * 2, points);

	for (int i = 1; i < leftSize - 1; i += 2)
		MergePoints(points, leftPos + i * offset, leftPos + (i + 1) * offset);

	if (leftSize % 2 == 0)
		MergePoints(points, leftPos + (leftSize - 1) * offset, rightPos);

	for (int i = 1 - leftSize % 2; i < rightSize - 1; i += 2)
		MergePoints(points, rightPos + i * offset, rightPos + (i + 1) * offset);
}

void ParallelBetcherSorter::Sort(int position, int length, vector<Point> &points) {
	if (length < 2)
		return;

	int middle = length / 2;
	Sort(position, middle, points);
	Sort(position + middle, length - middle, points);
	Merge(position, position + middle, middle, length - middle, 1, points);
}

bool ParallelBetcherSorter::Check(vector<Point> &points) {
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

void ParallelBetcherSorter::Sort(vector<Point> &points, SortAlgorithm algorithm) {
    if (algorithm == HeapSortAlgorithm) {
        HeapSort(points);
    }
    else if (algorithm == MergeSortAlgorithm) {
        MergeSort(points, 0, points.size() - 1);
    }
    else {
        sort(points.begin(), points.end(), comparator);
    }

    Sort(0, size, points);
}


vector<Point> GetPoints(int n1, int n2, int rank, int count, FillMode mode) {
    vector<Point> points(count);

    for (int index = 0; index < count; index++) {
        int i = (index + rank * count) / n2;
        int j = (index + rank * count) % n2;

        points[index] = GetPoint(i, j, n1, n2, mode);
    }

    return points;
}

void PrintPoints(const vector<Point> &points) {
    for (int i = 0; i < points.size(); i++)
        cout << points[i] << " ";
        
    cout << endl;
}

void Help() {
    cout << "Usage: ./task2 n1 n2 [coord] [mode] [path] [algorithm] [debug]" << endl;
    cout << "Arguments:" << endl;
    cout << "  n1        - first size of grid" << endl;
    cout << "  n2        - second size of grid" << endl;
    cout << "  coord     - coordinate for sorting (0 - x, 1 - y), default is 0 (x)" << endl;
    cout << "  mode      - mode of filling points, default is trigonometry" << endl;
    cout << "  path      - path to file for saving results, default 'results.jsonl'" << endl;
    cout << "  algorithm - sorting algorithm, default std::sort" << endl;
    cout << "  debug     - debug mode for printing points, not used by default" << endl << endl;

    cout << "Filling modes:" << endl;
    cout << "  a - ascending (x = y = index)" << endl;
    cout << "  d - descending (x = y = size - index - 1)" << endl;
    cout << "  r - random" << endl;
    cout << "  t - trigonometry (x = sin(phi), y = cos(phi), where phi = 2pi * index / size" << endl << endl;

    cout << "Algorithms:" << endl;
    cout << "  std       - standard sort (std::sort)" << endl;
    cout << "  mergesort - merge sort algorithm" << endl;
    cout << "  heapsort  - heap sort algorithm" << endl;
}

int main(int argc, const char** argv) {
    if (argc == 1 || (argc == 2 && (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")))) {
        Help();
        return 0;
    }

    int size, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(rank);

    int n1 = atoi(argv[1]);
    int n2 = atoi(argv[2]);
    int sortedCoord = argc <= 3 ? 0 : atoi(argv[3]);
    FillMode mode = argc <= 4 ? Trigonometry : GetFillMode(argv[4]);
    const char* path = argc <= 5 ? "results.jsonl" : argv[5];
    SortAlgorithm algorithm = argc <= 6 ? StdSortAlgorithm : GetAlgorithm(argv[6]);
    bool debug = argc == 8;

    int n = n1 * n2;
    int count = (n + size - 1) / size;

    vector<Point> points = GetPoints(n1, n2, rank, count, mode);

    if (debug) {
        cout << "Input points " << rank << ": ";
        PrintPoints(points);    
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    double start = MPI_Wtime();
    ParallelBetcherSorter sorter(rank, size, sortedCoord);
    sorter.Sort(points, algorithm);
    
    double delta = MPI_Wtime() - start;
    double time;
    MPI_Reduce(&delta, &time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (debug) {
        cout << "Sorted points " << rank << ": "; 
        PrintPoints(points);
    }

    bool isSorted = sorter.Check(points);

    if (rank == 0) {
        cout << "n1: " << n1 << endl;
        cout << "n2: " << n2 << endl;
        cout << "processors: " << size << endl;
        cout << "sorted coord: " << (sortedCoord == 0 ? "x" : "y") << endl;
        cout << "time: " << time << endl;
        cout << "status: " << (isSorted ? "sorted" : "not sorted") << endl;

        ofstream fout(path, ios::app);
        fout << "{\"p\": " << size << ", \"n1\": " << n1 << ", \"n2\": " << n2 << ", \"time\": " << time << "}" << endl;
        fout.close();
    }

    MPI_Finalize();
    return 0;
}