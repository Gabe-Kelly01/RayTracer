#ifndef SHAPE_HDR
#define SHAPE_HDR

#include "Tuple.h"
#include "Rgb.h"

class Shape {
    public:
        int type; // 1 indicates a plane, 2 is a sphere
        int materialType; // 1 indicates matte, 2 is metallic, 3 is transparent
        double k; // shininess factor if metallic, refractive index if transparent
        Tuple position;
        Rgb rDiff;
        Rgb rSpec;
        Rgb rAmb;

};

#endif