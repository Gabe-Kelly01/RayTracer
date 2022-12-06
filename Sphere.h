#include "Shape.h"
#include "Intersection.h"

class Sphere : public Shape {
    public:
        double radius;

        Sphere(Tuple pos, double rad, float reflectCoef, Rgb diff,
               Rgb spec, Rgb amb);

};