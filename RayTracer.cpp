// CS360 Project 5 - Gabe Kelly
#include <vector>
#include <limits>
#include <algorithm>

#include "Lighting.h"
#include "Intersection.h"
#include "Plane.h"
#include "Sphere.h"
#include "easyppm.h"

#define MAX_RECURSION_DEPTH 5

struct Ray {
    Tuple point;
    Tuple direction;
};

struct Node {
    Shape* shapePtr;
    const int shapeType;
};

struct LightSrc {
    Tuple position;
    Rgb iAmb;
    Rgb iDiff;
    Rgb iSpec;
    int specExp;
};

float calculateFresnel(const float &index, const Tuple &normalHit, const double &k) {
    return 0.0;
}

bool rayHitsPlane(const Tuple& rayOriginPoint, const Tuple& rayDirectionVector, const Plane& plane, Intersection& intersect) {
    double t;
    double rayDotNorm = rayDirectionVector.dot(plane.normVector);

    if (rayDotNorm != 0.0) {
        t = (plane.position - rayOriginPoint).dot(plane.normVector) / rayDotNorm;
        if (t < 0) {
            return false;
        }

        intersect = Intersection((Shape *) &plane, t);
        return true;
    }

    return false;
}

bool rayHitsSphere( const Tuple& rayOriginPoint, const Tuple& rayDirectionVector, const Sphere& sphere, Intersection& intersect) {
    // Simplify variables for quadratic computation
    Tuple V = rayOriginPoint - sphere.position;

    double A = rayDirectionVector.dot(rayDirectionVector);
    double B = 2 * V.dot(rayDirectionVector);
    double C = V.dot(V) - pow(sphere.radius, 2);
    double discriminant = pow(B, 2) - (4 * A * C);

    if (discriminant < 0) {
        return false;
    } else {
        double t_1 = (-B - sqrt(discriminant)) / 2 * A;
        double t_2 = (-B + sqrt(discriminant)) / 2 * A;
        double t;


        if (t_1 < 0) {
            t = t_2;

            if (t_2 < 0) {
                return false;
            }
        }

        t = t_1;

        intersect = Intersection((Shape *) &sphere, t);
        return true;
    }
}

bool inShadow(Tuple& intersectPoint, Tuple& lightPoint, const std::vector<Node>& scene) {
    intersectPoint = intersectPoint + (0.001 * intersectPoint - lightPoint);
    Tuple intersectToLight = lightPoint - intersectPoint;
    double dL = intersectToLight.magnitude();
    intersectToLight.normalize();

    Intersection it;
    double t;
    for (Node n : scene) {
        if (n.shapeType == 1) {
            auto obj = (Plane *) n.shapePtr;

            if (rayHitsPlane(intersectPoint, intersectToLight, *obj, it)) {
                t = it.tHit;
                
                if (t < dL) {
                    return true;
                }
            }
        } else {
            auto obj = (Sphere *) n.shapePtr;

            if (rayHitsSphere(intersectPoint, intersectToLight, *obj, it)) {
                t = it.tHit;
                
                if (t < dL) {
                    return true;
                }
            }
        }
    }

    return false;
}

Rgb Trace(Ray& ray, const std::vector<Node>& scene, const std::vector<LightSrc>& lSource, int recursion_depth) {
    // Declare intersection vector
    std::vector<Intersection> hits;

    // Iterate over every object in scene
    for (Node n : scene) {
        Intersection it;

        // Cast node to correct shape
        // shapeType == 1 -> plane
        // shapeType == 2 -> sphere
        if (n.shapeType == 1) {
            auto obj = (Plane *) n.shapePtr;

            if (rayHitsPlane(ray.point, ray.direction, *obj, it)) {
                hits.push_back(it);
            }
        } else {
            auto obj = (Sphere *) n.shapePtr;

            if (rayHitsSphere(ray.point, ray.direction, *obj, it)) {
                hits.push_back(it);
            }
        }
    }

    // Find the closest intersection, unpack shape and find point of intersection
    sort(hits.begin(), hits.end(), sortByHT);

    if (hits.size() == 0) {
        return Rgb();
    }
    Intersection minHit = hits[0];
    Shape curObj = *minHit.obj;
    Tuple intersectPoint = ray.point + ray.direction * minHit.tHit;

    // Calculate lighting from all light sources in the scene
    Tuple lightNormal;
    if (curObj.type == 1) {
        lightNormal = Tuple(0, 0, -1);

    } else {
        lightNormal = intersectPoint - curObj.position;
        lightNormal.normalize();

        if (ray.direction.dot(lightNormal) < 0) {
            lightNormal.z *= -1;
        }
    }

    // Rgb totalLight = Rgb();
    // // Calculate for a reflective or transparent object
    // if ((curObj.reflectiveness != 1) and recursion_depth < MAX_RECURSION_DEPTH) {
    //     // Compute reflected ray
    //     Tuple reflectDirection = ray.direction - (lightNormal * 2 * ray.direction.dot(lightNormal));
    //     reflectDirection.normalize();

    //     Tuple reflectOrigin = intersectPoint + (lightNormal * 0.001);
    //     Ray reflectionRay = {reflectOrigin, reflectDirection};
    //     totalLight = totalLight + (curObj.reflectiveness * Trace(reflectionRay, scene, lSource, recursion_depth + 1));
    // } else {
    //     // Calculate for a matte object
    //     for (LightSrc src : lSource) {
    //         // Compute ambient light regardless, only compute specular and diffuse lighting if point is not in shadows
    //         totalLight = totalLight + lightAmbient(curObj.rAmb, src.iAmb);
    //         bool hasShadow = inShadow(intersectPoint, src.position, scene);
    //         if (not hasShadow) {
    //             totalLight = totalLight +
    //                          lightDiffuse(curObj.rDiff, curObj.position, lightNormal, src.iDiff, src.position);
    //             totalLight = totalLight +
    //                          lightSpecular(curObj.rSpec, curObj.position, lightNormal, src.iSpec, src.position,
    //                                        ray.point, src.specExp);
    //         }
    //     }
    // }

    Rgb totalLight = Rgb();
    for (LightSrc src : lSource) {

        if ((curObj.reflectiveness != 1) and (recursion_depth < MAX_RECURSION_DEPTH)) {
            // Compute reflected ray
            Tuple reflectDirection = ray.direction + (lightNormal * 2 * ray.direction.dot(lightNormal));
            reflectDirection.normalize();

            Tuple reflectOrigin = intersectPoint + (lightNormal * 0.001);
            Ray reflectionRay = {reflectOrigin, reflectDirection};
            totalLight = totalLight + (curObj.reflectiveness * Trace(reflectionRay, scene, lSource, recursion_depth + 1));
        }
        else {
            totalLight = totalLight + lightAmbient(curObj.rAmb, src.iAmb);
            
            bool hasShadow = inShadow(intersectPoint, src.position, scene);
            if (not hasShadow) {
                totalLight = totalLight +
                            lightDiffuse(curObj.rDiff, curObj.position, lightNormal, src.iDiff, src.position);
                totalLight = totalLight +
                            lightSpecular(curObj.rSpec, curObj.position, lightNormal, src.iSpec, src.position,
                                            ray.point, src.specExp);
            }
        }

        // if ((curObj.reflectiveness != 1) and (recursion_depth < MAX_RECURSION_DEPTH)) {
        //     // Compute reflected ray
        //     Tuple reflectDirection = ray.direction + (lightNormal * 2 * ray.direction.dot(lightNormal));
        //     reflectDirection.normalize();

        //     Tuple reflectOrigin = intersectPoint + (lightNormal * 0.001);
        //     Ray reflectionRay = {reflectOrigin, reflectDirection};
        //     totalLight = totalLight + (curObj.reflectiveness * Trace(reflectionRay, scene, lSource, recursion_depth + 1));
        // }
    }

    return totalLight;
}
#pragma clang diagnostic pop

void Render(const std::vector<Node>& scene, const std::vector<LightSrc>& lSource) {
    double front_clip = 6.0;
    double w = 6.0;
    double h = 6.0;

    // Bottom-left point of image output window
    Tuple b_left = Tuple(-w/2, -h/2, front_clip, 1);

    // Point/Vector constants that will be reused in the main loop
    Tuple X = Tuple(1.0, 0, 0);
    Tuple Y = Tuple(0, 1.0, 0);
    Tuple cameraPoint = Tuple(0,0,0, 1);

    int image_pixel_size = 750;
    PPM ray_trace_image = easyppm_create(image_pixel_size, image_pixel_size, IMAGETYPE_PPM);
    easyppm_clear(&ray_trace_image, easyppm_rgb(255, 255, 255));

    for (double t = 0.0; t < h; t += (h / image_pixel_size)) {
        for (double s = 0.0; s < w; s += (w / image_pixel_size)) {
            // Calculate point P
            Tuple P = b_left + s * X + t * Y;

            // Create normalized ray from origin to P
            Tuple direction = P - cameraPoint;
            direction.normalize();
            Ray ray = {cameraPoint, direction};
            Rgb totalLight = Trace(ray, scene, lSource, 0);

            // Shade pixel
            int x = (int)((s * image_pixel_size / w) * (1.0 + std::numeric_limits<double>::epsilon()));
            int y = (int)((t * image_pixel_size / h) * (1.0 + std::numeric_limits<double>::epsilon()));
            easyppm_set(&ray_trace_image, x, y, easyppm_rgb(255 * totalLight.getR(),
                                                                    255 * totalLight.getG(),
                                                                    255 * totalLight.getB()));
        }
    }

    // Write PPM image to file
    easyppm_write(&ray_trace_image, "image.ppm");
    easyppm_destroy(&ray_trace_image);
}


int main() {
    Plane p1 = Plane(Tuple(0, -4, 15, 1), Tuple(0, 0, -1), Rgb(1, 0, 0),
                     Rgb(1, 0, 0), Rgb(1, 0, 0));
    Node o1 = {&p1, p1.type};

    Plane p2 = Plane(Tuple(0, -5, 15, 1), Tuple(0, 1, 0), Rgb(0, 0, 1),
                     Rgb(0, 0, 0.5), Rgb(0, 0, 1));
    Node o2 = {&p2, p2.type};

    Sphere s1 = Sphere(Tuple(1.5, 1.5, 6, 1), 1, 1, Rgb(0, 1, 0),
                       Rgb(0, 1, 0), Rgb(0, 1, 0));
    Node o3 = {&s1, s1.type};

    Sphere s2 = Sphere(Tuple(0, 1, 10, 1), 2, 0.4, Rgb(0.8,0.6,0.8),
                       Rgb(0.8,0.6,0.8), Rgb(0.8,0.6,0.8));
    Node o4 = {&s2, s2.type};

    Sphere s3 = Sphere(Tuple(-1.5, -1.5, 6, 1), 1, 1, Rgb(0.5,0,0.5),
                       Rgb(0.5,0,0.5), Rgb(0.5,0,0.5));
    Node o5 = {&s3, s3.type};

    Sphere s4 = Sphere(Tuple(-1.5, 1.5, 6, 1), 1, 1, Rgb(0,0.5,0.5),
                       Rgb(0,0.5,0.5), Rgb(0,0.5,0.5));
    Node o6 = {&s4, s4.type};

    Sphere s5 = Sphere(Tuple(2.0, -1.5, 14, 1), 2, 1, Rgb(0.2,0.6,0.8),
                       Rgb(0.2,0.6,0.8), Rgb(0.2,0.6,0.8));
    Node o7 = {&s5, s5.type};

    std::vector<Node> scene {o1, o2, o3, o4, o5, o6, o7};

    LightSrc l1 = {Tuple(3, 3, 6, 1),
                   Rgb(0.3, 0.3, 0.3),
                   Rgb(0.5, 0.5, 0.5),
                   Rgb(0.8, 0.8, 0.8), 2};
    LightSrc l2 = {Tuple(-3, 3, 6, 1),
                   Rgb(0.3, 0.3, 0.3),
                   Rgb(0.2, 0.2, 0.2),
                   Rgb(0.4, 0.4, 0.4), 1};

    std::vector<LightSrc> lights {l1, l2};

    Render(scene, lights);

    return 0;
}

void GenerateWaveAnimation() {
    LightSrc l1 = {Tuple(0, -3, 6, 1),
                   Rgb(0.3, 0.3, 0.3),
                   Rgb(0.5, 0.5, 0.5),
                   Rgb(0.8, 0.8, 0.8), 2};
    std::vector<LightSrc> lights {l1};

    const double dMax = 3;
    int numSpheres = 5;
    int frames = 60;

    double t = 0.0;
    for (int fc = 0; fc <= frames; fc++) {
        std::vector<Sphere> spheres;
        for (int s = 0; s < numSpheres; s++) {
            void();
        }
    }

}

double getHeight(double t, double dMax) {
    return dMax - (dMax * sin(t));
}
