#pragma once
#include <iostream>


enum SortAlgorithm {
    MergeSortAlgorithm,
    HeapSortAlgorithm,
    StdSortAlgorithm
};

struct Point {
    float coord[2];
    int index;
};


SortAlgorithm GetAlgorithm(const char *arg) {
    if (!strcmp(arg, "heapsort"))
        return HeapSortAlgorithm;

    if (!strcmp(arg, "mergesort"))
        return MergeSortAlgorithm;

    if (!strcmp(arg, "std"))
        return StdSortAlgorithm;

    std::cout << "Warning: invalid sorting algorithm \"" << arg << "\". Used std" << std::endl;

    return StdSortAlgorithm;
}

std::ostream& operator <<(std::ostream& os, const Point& point) {
    return os << point.index << ": {" << point.coord[0] << " " << point.coord[1] << "}";
}
