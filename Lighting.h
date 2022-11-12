#ifndef LIGHTING_HDR
#define LIGHTING_HDR

#include "Rgb.h"
#include "Tuple.h"

/*
 * lightAmbient
 *
 * Compute and return ambient reflectance as an Rgb value, where 0.0 < r, g, b <= 1.0
 * materialReflectance - Given material RGB reflectance property - what fraction of incoming red, green, blue light
 *                       is reflected by this material. Type Rgb contains three type double values each between 0.0 and 1.0.
 * lightIntensity - intensity of ambient light given as RGB coefficients. Yes, lights can have color.
 *                       A maximally bright white light will have Rgb coefficients of 1.0, 1.0, 1.0.
 * Return computed amount of ambient light reflected by given material as an Rgb.
 *        Intensity of reflected light = materialReflectance * lightIntensity
 */
Rgb lightAmbient( const Rgb& materialReflectance, const Rgb& lightIntensity );

/*
 * lightDiffuse
 *
 * Compute and return diffuse reflectance as an Rgb value, where 0.0 < r, g, b <= 1.0
 * materialReflectance - Given material RGB reflectance property - what fraction of incoming red, green, blue light
 *                       is reflected by this material. Type Rgb contains three type double values each between 0.0 and 1.0.
 * objectPoint - position of point on object that we wish to shade
 * objectNormal - unit-length normal vector perpendicular to surface of material at point P
 * lightIntensity - intensity of ambient light given as RGB coefficients. Yes, lights can have color.
 *                  A maximally bright white light will have Rgb coefficients of 1.0, 1.0, 1.0.
 * lightPoint - position of a light source.
 * Return computed amount of ambient light reflected by given material as an Rgb.
 *        Intensity of reflected light = materialReflectance * lightIntensity * ( L DOT N ),
 *        where L is the unit-length vector from objectPoint to lightPoint and
 *              N is the unit-length normal vector that is perpendicular to the material surface at objectPoint.
 */
Rgb lightDiffuse(const Rgb& materialReflectance, const Tuple& objectPoint, const Tuple& objectNormal, const Rgb& lightIntensity, const Tuple& lightPoint);

/*
 * lightSpecular
 *
 * Compute and return specular reflectance as an Rgb value, where 0.0 < r, g, b <= 1.0
 * materialReflectance - Given material RGB reflectance property - what fraction of incoming red, green, blue light
 *                       is reflected by this material. Type Rgb contains three type double values each between 0.0 and 1.0.
 * objectPoint - position of point on object that we wish to shade
 * objectNormal - unit-length normal vector perpendicular to surface of material at point P
 * lightIntensity - intensity of ambient light given as RGB coefficients. Yes, lights can have color.
 *                  A maximally bright white light will have Rgb coefficients of 1.0, 1.0, 1.0.
 * lightPoint - position of a light source.
 * eyePoint - position of a eye or virtual camera
 * exponent - specular exponent given as a positive integer value >= 1.
 * Return computed amount of specular light reflected by given material as an Rgb.
 *        Intensity of reflected light = materialReflectance * lightIntensity * ( RL DOT EV )^exponent,
 *        where RL is the unit-length vector that reflects the vector from objectPoint to lightPoint about objectNormal vector
 *              EV is the unit-length normal vector from the objectPoint to the eye point,
 *              ^exponent means to raise the quantity to the exponent power. Note that C++ does not support the ^ operator.
 */
Rgb lightSpecular(const Rgb& materialReflectance, const Tuple& objectPoint, const Tuple& objectNormal, 
                  const Rgb& lightIntensity, const Tuple& lightPoint,
									const Tuple& eyePoint, int exponent);

#endif