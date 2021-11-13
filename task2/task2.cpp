#include <iostream>
#include <vector>
#include <mpi.h>

#include "Entities.hpp"
#include "PointGenertor.hpp"

using namespace std;


vector<Point> GetPoints(int n1, int n2, int rank, int count, FillMode mode) {
    vector<Point> array(count);

    for (int index = 0; index < count; index++) {
        int i = (index + rank * count) / n2;
        int j = (index + rank * count) % n2;

        array[index] = GetPoint(i, j, n1, n2, mode);
    }

    return array;
}


int main(int argc, const char** argv) {
    int size, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(rank);

    int n1 = atoi(argv[1]);
    int n2 = atoi(argv[2]);
    int sortedCoord = argc <= 3 ? 0 : atoi(argv[3]);
    FillMode mode = argc <= 4 ? Trigonometry : GetFillMode(argv[4]);

    int n = n1 * n2;
    int count = (n + size - 1) / size;

    vector<Point> array = GetPoints(n1, n2, rank, count, mode);
    cout << "Input array " << rank << ": ";

    for (int i = 0; i < count; i++) {
        cout << array[i] << " ";
    }

    cout << endl;

    // cout << endl << "Sorted array: ";
    // for (int i = 0; i < n; i++)
    //     cout << array[i] << " ";
    MPI_Finalize();
    return 0;
}