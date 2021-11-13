#include <iostream>
#include <vector>
#include <algorithm>

#include <mpi.h>

#include "Entities.hpp"
#include "PointGenertor.hpp"

using namespace std;

int sortedCoord; // потому что на bluegene нельзя в лямбды ((((

MPI_Datatype MakePointType() {
    Point point;
    MPI_Datatype pointType;
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

    return pointType;
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

bool ComparePoints(Point a, Point b) {
    return a.coord[sortedCoord] < b.coord[sortedCoord];
}

void MakeMergePoints(const vector<Point> &pointsI, const vector<Point> &pointsJ, vector<Point> &points) {
    size_t i = 0;
    size_t j = 0;
    size_t k = 0;

    while (i < pointsI.size() && j < pointsJ.size())
        points[k++] = ComparePoints(pointsI[i], pointsJ[j]) ? pointsI[i++] : pointsJ[j++];

    while (i < pointsI.size())
        points[k++] = pointsI[i++];

    while (j < pointsJ.size())
        points[k++] = pointsJ[j++];
}

void MergePoints(vector<Point> &pointsI, int rank, int pi, int pj, MPI_Datatype pointType) {
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

    cout << rank << ": merge with " << pj << endl; 
}

void Merge(int leftPos, int rightPos, int leftSize, int rightSize, int offset, vector<Point> &points, int rank, MPI_Datatype pointType) {
	if (leftSize == 0 || rightSize == 0)
		return;

	if (leftSize == 1 && rightSize == 1) {
		MergePoints(points, rank, leftPos, rightPos, pointType);
		return;
	}

	Merge(leftPos, rightPos, (leftSize + 1) / 2, (rightSize + 1) / 2, offset * 2, points, rank, pointType);
	Merge(leftPos + offset, rightPos + offset, leftSize / 2, rightSize / 2, offset * 2, points, rank, pointType);

	for (int i = 1; i < leftSize - 1; i += 2)
		MergePoints(points, rank, leftPos + i * offset, leftPos + (i + 1) * offset, pointType);

	if (leftSize % 2 == 0)
		MergePoints(points, rank, leftPos + (leftSize - 1) * offset, rightPos, pointType);

	for (int i = 1 - leftSize % 2; i < rightSize - 1; i += 2)
		MergePoints(points, rank, rightPos + i * offset, rightPos + (i + 1) * offset, pointType);
}

void Sort(int position, int size, vector<Point> &points, int rank, MPI_Datatype pointType) {
	if (size < 2)
		return;

	int middle = size / 2;
	Sort(position, middle, points, rank, pointType);
	Sort(position + middle, size - middle, points, rank, pointType);
	Merge(position, position + middle, middle, size - middle, 1, points, rank, pointType);
}


int main(int argc, const char** argv) {
    int size, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Datatype pointType = MakePointType();

    srand(rank);

    int n1 = atoi(argv[1]);
    int n2 = atoi(argv[2]);
    sortedCoord = argc <= 3 ? 0 : atoi(argv[3]);
    FillMode mode = argc <= 4 ? Trigonometry : GetFillMode(argv[4]);
    FillMode mode = argc <= 4 ? Trigonometry : GetFillMode(argv[4]);

    int n = n1 * n2;
    int count = (n + size - 1) / size;

    vector<Point> points = GetPoints(n1, n2, rank, count, mode);
    cout << "Input points " << rank << ": ";
    PrintPoints(points);

    sort(points.begin(), points.end(), ComparePoints);

    cout << "Sorted points " << rank << ": ";
    PrintPoints(points);

    Sort(0, size, points, rank, pointType);
    
    cout << "Merged points " << rank << ": "; 
    PrintPoints(points);

    MPI_Finalize();
    return 0;
}