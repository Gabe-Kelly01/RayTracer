#ifndef INTERSECTION_HDR
#define INTERSECTION_HDR

#include "Shape.h"

// Class representing intersections
class Intersection{
public:
    Shape* obj;
    double tHit;

    Intersection();
    Intersection(Shape* obj, double hit);
};

bool sortByHT(const Intersection& a, const Intersection& b);

#endif