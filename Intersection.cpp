#include "Intersection.h"

Intersection::Intersection() {
    this->obj = nullptr;
    this->tHit = 0.0;
}

Intersection::Intersection(Shape *obj, double hit) {
    this->obj = obj;
    this->tHit = hit;
}

bool sortByHT(const Intersection& a, const Intersection& b) {
    return a.tHit <= b.tHit;
}