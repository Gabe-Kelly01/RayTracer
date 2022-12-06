#include "Sphere.h"

Sphere::Sphere(Tuple pos, double rad, float reflectCoef, Rgb diff, Rgb spec, Rgb amb) {
    this->type = 2;
    this->position = pos;
    this->reflectiveness = reflectCoef;
    this->radius = rad;
    this->rDiff = diff;
    this->rSpec = spec;
    this->rAmb = amb;
}