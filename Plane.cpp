#include "Plane.h"

Plane::Plane(Tuple pos, Tuple normVector, Rgb diff, Rgb spec, Rgb amb) {
    this->type = 1;
    this->position = pos;
    this->normVector = normVector;
    this->rDiff = diff;
    this->rSpec = spec;
    this->rAmb = amb;
}