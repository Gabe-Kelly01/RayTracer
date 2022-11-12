// CS360 Project 4 - Gabe Kelly
#include "Lighting.h"

Rgb lightAmbient(const Rgb& materialReflectance, const Rgb& lightIntensity) {
    return materialReflectance * lightIntensity;
}

Rgb lightDiffuse(const Rgb& materialReflectance, const Tuple& objectPoint, const Tuple& objectNormal, const Rgb& lightIntensity, const Tuple& lightPoint) {
    Tuple L = lightPoint - objectPoint;
    L.normalize();
    double incidenceAngle = L.dot(objectNormal);

    return materialReflectance * lightIntensity * incidenceAngle;
}

Rgb lightSpecular(const Rgb& materialReflectance, const Tuple& objectPoint, const Tuple& objectNormal,
                  const Rgb& lightIntensity, const Tuple& lightPoint,
                  const Tuple& eyePoint, int exponent) {
    Tuple L = objectPoint - lightPoint;
    L.normalize();

    Tuple RV = 2 * L.dot(objectNormal) * objectNormal - L;
    RV.normalize();

    Tuple EV = eyePoint - objectPoint;
    EV.normalize();

    double prod = EV.dot(RV);
    double alpha = pow(prod, exponent);

//    std::cout << "L is: " << L << std::endl;
//    std::cout << "RV is: " << RV << std::endl;
//    std::cout << "EV is: " << EV << std::endl;
//    std::cout << "EV dot RV is: " << prod << std::endl;
//    std::cout << "alpha is: " << alpha << std::endl;
    return materialReflectance * lightIntensity * alpha;

}