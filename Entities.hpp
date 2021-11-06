#pragma once
#include <iostream>


struct Point {
    float coord[2];
    int index;
};

std::ostream& operator <<(std::ostream& os, const Point& point) {
    return os << point.index << ": {" << point.coord[0] << " " << point.coord[1] << "}";
}
