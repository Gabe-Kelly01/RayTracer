#include "Sphere.h"

Sphere::Sphere(Tuple pos, double rad, Rgb diff, Rgb spec, Rgb amb) {
    this->type = 2;
    this->position = pos;
    this->radius = rad;
    this->rDiff = diff;
    this->rSpec = spec;
    this->rAmb = amb;
}