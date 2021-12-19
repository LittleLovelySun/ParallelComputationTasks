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


int GetCuttedEdges(vector<Point> &points, int n1, int n2) {
    int n = n1 * n2;
    int edges = 0;

    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (Distance(points[i], points[j], n2) == 1 && points[i].domain != points[j].domain)
                edges++;

    return edges;
}


// TODO
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
    GrigType type = argc <= 4 ? Grid : GetGridType(argv[4]);
    const char* path = argc <= 5 ? NULL : argv[5];

    GridGenerator generator(size, rank, type);
    vector<Point> points = generator.Generate(n1, n2);

    MPI_Barrier(MPI_COMM_WORLD);

    double start = MPI_Wtime();
    GridCoordinateDecomposer decomposer(size, rank);
    decomposer.Decompose(points, k);

    double delta = MPI_Wtime() - start;
    double time;
    MPI_Reduce(&delta, &time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    vector<Point> totalPoints = decomposer.JoinPoints(points);

    if (rank == 0) {
        if (path) {
            ofstream f(path);
            PrintPoints(totalPoints, n2, f);
            f.close();
        }

        int edges = GetCuttedEdges(totalPoints, n1, n2);

        cout << "n1: " << n1 << endl;
        cout << "n2: " << n2 << endl;
        cout << "k: " << k << endl;
        cout << "processors: " << size << endl;
        cout << "edges: " << edges << endl;
        cout << "time: " << time << endl;
    }

    MPI_Finalize();
    return 0;
}