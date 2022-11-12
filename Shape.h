#ifndef SHAPE_HDR
#define SHAPE_HDR

#include "Tuple.h"
#include "Rgb.h"

class Shape {
    public:
        int type; // 1 indicates a plane, 2 is a sphere
        Tuple position;
        Rgb rDiff;
        Rgb rSpec;
        Rgb rAmb;

};

#endif