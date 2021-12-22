#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <cstring>

#include <mpi.h>

#include "Entities.hpp"
#include "GridGenerator.hpp"
#include "GridCoordinateDecomposer.hpp"

using namespace std;


void Help() {
    cout << "Usage: ./task2 n1 n2 k [grid type] [decompose type] [path]" << endl;
    cout << "Arguments:" << endl;
    cout << "  n1              - first size of grid" << endl;
    cout << "  n2              - second size of grid" << endl;
    cout << "  k               - number of domains for decomposing" << endl;
    cout << "  grid type       - type of grid for generate points, default is 'grid'" << endl;
    cout << "  decompose type  - sorting algorithm, default std::sort" << endl;
    cout << "  path            - path to file for saving results" << endl;

    cout << "Grid types:" << endl;
    cout << "  g - rectangular grid (x = i, y = j)" << endl;
    cout << "  s - spiral" << endl;
    cout << "  c - circle" << endl;
    cout << "  f - flower" << endl << endl;

    cout << "Decompose types:" << endl;
    cout << "  c - sequential change axis" << endl;
    cout << "  i - minimizing intersections on every step" << endl;
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

    int n1 = atoi(argv[1]);
    int n2 = atoi(argv[2]);
    int k = atoi(argv[3]);
    GridType gridType = argc <= 4 ? Grid : GetGridType(argv[4]);
    DecomposeType decomposeType = argc <= 5 ? Intersections : GetDecomposeType(argv[5]);
    const char* path = argc <= 6 ? NULL : argv[6];

    GridGenerator generator(size, rank, gridType);
    vector<Point> points = generator.Generate(n1, n2);

    MPI_Barrier(MPI_COMM_WORLD);

    double start = MPI_Wtime();
    GridCoordinateDecomposer decomposer(size, rank, n1, n2);
    decomposer.Decompose(points, k, decomposeType);

    double delta = MPI_Wtime() - start;
    double time;
    MPI_Reduce(&delta, &time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (path)
        decomposer.Save(points, path);

    int edges = decomposer.GetTotalEdges(points, decomposeType);

    if (rank == 0) {
        ofstream fout("result.jsonl", ios::app);
        fout << "{";
        fout << "\"n1\": " << n1;
        fout << ", \"n2\": " << n2;
        fout << ", \"grid\": \"" << gridType << "\"";
        fout << ", \"method\": \"" << decomposeType << "\"";
        fout << ", \"k\": " << k;
        fout << ", \"p\": " << size;
        fout << ", \"time\": " << time;
        fout << ", \"edges\": " << edges;
        fout << "}" << endl;
        fout.close();
    }

    MPI_Finalize();
    return 0;
}