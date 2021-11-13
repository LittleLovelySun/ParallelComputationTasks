#include <iostream>
#include <vector>

#include "Entities.hpp"
#include "PointGenertor.hpp"

using namespace std;


int main(int argc, const char** argv) {
    int n1 = atoi(argv[1]);
    int n2 = atoi(argv[2]);
    int sortedCoord = argc <= 3 ? 0 : atoi(argv[3]);
    FillMode mode = argc <= 4 ? Random : GetFillMode(argv[4]);

    int n = n1 * n2;

    vector<Point> array(n);
    cout << "Input array  : ";

    for (int i = 0; i < n; i++) {
        array[i] = GetPoint(i / n2, i % n2, n1, n2, mode);
        cout << array[i] << " ";
    }

    cout << endl << "Sorted array: ";
    for (int i = 0; i < n; i++)
        cout << array[i] << " ";

    return 0;
}