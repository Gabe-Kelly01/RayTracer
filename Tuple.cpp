#include "Tuple.h"

// Tuple has format (x,y,z,w) where x, y, and z represent the coordinate/vector data
// w distinguishes whether a tuple represents a vector or a point, I have changed it to an
// int for less of a headache
// w = 0 => tuple is a vector
// w = 1 => tuple is a point

// Default constructor for when arguments are not provided
// Initializes a new point tuple with value (0.0, 0.0, 0.0, 0)
Tuple::Tuple() {
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
    this->w = 0;
}

// Constructor that creates a Tuple representing a vector (x,y,z,0)
Tuple::Tuple(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = 0;
}

// Constructor that creates a Tuple representing a point (x,y,z,1)
// w is expected to be a positive value but will be stored as 1 regardless
Tuple::Tuple(double x, double y, double z, int w) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = 1;
}

// Copy constructor that initializes a new Tuple with identical value to the source Tuple
Tuple::Tuple(const Tuple& source) {
    this->x = source.x;
    this->y = source.y;
    this->z = source.z;
    this->w = source.w;
}

// Returns whether this Tuple is vector
// Abs method is not used because I opted to use an int as a flag
bool Tuple::isVector() const {
    return this->w == 0;
}

// Returns whether this Tuple is a point
// Simply the negation of isVector()
bool Tuple::isPoint() const {
    return this->w == 1;
}

// Increments this Tuple's values by the values of otherTuple as follows
// (x,y,z) + (a,b,c) = (x+a, y+b, z+c)
// The add method adheres to the following conversion rules:
// point + vector = point
// point + point = point
// vector + vector = vector
void Tuple::add(const Tuple &otherTuple) {
    this->x += otherTuple.x;
    this->y += otherTuple.y;
    this->z += otherTuple.z;

    // Check if point-to-vector conversion is needed
    if (this->isPoint() and otherTuple.isVector() or this->isVector() and otherTuple.isPoint() ) {
        this->w = 1;
    }
}

// Decrements this tuples values by the value of otherTuple as follows:
// (x,y,z) - (a,b,c) = (x-a,y-b,z-c)
// The sub method adheres to the following conversion rules:
// point - vector = point
// point - point = vector
// vector - vector = vector
void Tuple::sub(const Tuple &otherTuple) {
    this->x -= otherTuple.x;
    this->y -= otherTuple.y;
    this->z -= otherTuple.z;

    // Handling vector-point conversion
    if (this->isPoint() and otherTuple.isVector() or this->isVector() and otherTuple.isPoint()) {
        this->w = 1;
    } else if (this->isPoint() and otherTuple.isPoint()) {
        this->w = 0;
    }
}

// Multiplies the (x,y,z) value of the tuple by the scalar input S
// Result will be: (x,y,z,w) * S = (x*S, y*S, z*S, w)
void Tuple::multScalar(double S) {
    this->x *= S;
    this->y *= S;
    this->z *= S;
}

// Returns the magnitude (length) of the tuple
// If the tuple is a vector, this is the length of the vector
// If the tuple is a point, this is the distance of the point from (0,0,0)
// The w value is not affected
double Tuple::magnitude() const {
    // Taking powers of components separately for readability
    double x_r = pow(this->x, 2);
    double y_r = pow(this->y, 2);
    double z_r = pow(this->z, 2);

    // Using Pythagorean theorem to find the length
    // magnitude((x,y,z)) = sqrt(x^2 + y^2 + z^2)
    double result = sqrt(x_r + y_r + z_r);

    return result;
}

// If the tuple is vector, convert it to its unit vector, or a vector of
// length 1 that still points in the same direction. If the tuple is a point
// or if it has magnitude <0, then it will do nothing.
// Value of w is not changed
void Tuple::normalize() {
    if (this->isPoint()) {
        std::cout << "Cannot normalize a point (w=1)" << std::endl;
    } else {
        double mag = this->magnitude();

        if (mag <= 0) {
            std::cout << "Cannot normalize a vector with magnitude = 0" << std::endl;
        } else {
            this->x /= mag;
            this->y /= mag;
            this->z /= mag;
        }
    }
}

// Return dot product (a scalar) of this and otherTuple, only if they are both vectors
// If both tuples are not vectors, will return 0.0;
// Value of w is not changed
double Tuple::dot(const Tuple &otherTuple) const {
    if (this->isVector() and otherTuple.isVector()) {
        // Dot product of vectors (x,y,z) and (a,b,c)
        // = x*a + y*b + z*c
        double x_comp = this->x * otherTuple.x;
        double y_comp = this->y * otherTuple.y;
        double z_comp = this->z * otherTuple.z;
        double dot_prod = x_comp + y_comp + z_comp;

        return dot_prod;
    } else {
        return 0.0;
    }
}

void Tuple::print() const {
    std::cout << "(" << this->x << ", " << this->y << ", " << this->z << ", " << this->w << ")" << std::endl;
}

// Overload for the << operator to allow for Tuple objects to be printed like:
// Tuple A = Tuple(1.0, 2.0, 3.0, 1)
// std::cout << A << std::endl
// Result: (1.0, 2.0, 3.0, 1)
std::ostream& operator<<(std::ostream& os, const Tuple& T) {
    os << "(" << T.x << ", " << T.y << ", " << T.z << ", " << T.w << ")";
    return os;
}

// Overload for the addition operator. Allows for the addition of Tuple objects
// as you would with any other type. Creates copies of the input Tuples and returns
// new Tuple with the same value you would get with a.add(b)
Tuple operator+(const Tuple& a, const Tuple& b) {
    Tuple a_copy = Tuple(a);
    Tuple b_copy = Tuple(b);
    a_copy.add(b_copy);

    return a_copy;
}

// Overload for the subtraction operator. Allows for the subtraction of Tuple objects
// as you would with any other type. Creates copies of the input Tuples and returns
// a new Tuple with same value you would get with a.sub(b)
Tuple operator-(const Tuple& a, const Tuple& b) {
    Tuple a_copy = Tuple(a);
    Tuple b_copy = Tuple(b);
    a_copy.sub(b_copy);

    return a_copy;
}

// Overloads for the multiplication operator. Allows for scalar multiplication of
// Tuple objects using *, as you would with other types. Creates a copy of the input
// Tuple and returns a Tuple with the value a.multScalar(s). Two versions are provided
// so that both "a * s" and "s * a" will work in code.
Tuple operator*(const Tuple& a, double s) {
    Tuple copy = Tuple(a);
    copy.multScalar(s);

    return copy;
}

Tuple operator*(double s, const Tuple& a) {
    Tuple copy = Tuple(a);
    copy.multScalar(s);

    return copy;
}
