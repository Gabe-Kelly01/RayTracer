#include "Shape.h"

class Plane : public Shape {
    public:
        Tuple normVector;

        Plane(Tuple pos, Tuple normVector, Rgb diff,
              Rgb spec, Rgb amb);
};