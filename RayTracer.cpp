// CS360 Project 5 - Gabe Kelly
#include <vector>
#include <limits>
#include <algorithm>

#include "Lighting.h"
#include "Intersection.h"
#include "Plane.h"
#include "Sphere.h"
#include "easyppm.h"


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

bool rayHitsPlane(const Tuple& rayOriginPoint, const Tuple& rayDirectionVector, const Plane& plane, Intersection& intersect) {
    double t;
    double rayDotNorm = rayDirectionVector.dot(plane.normVector);

    if (rayDotNorm != 0.0) {
        t = (plane.position - rayOriginPoint).dot(plane.normVector) / rayDotNorm;
        if (t < 0) {
            return false;
        }
        std::cout << t << std::endl;
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


        if (t_1 < 0 && t_2 < 0) {
            return false;
        } else if ((t_1 > 0) != (t_2 > 0)) {
            t = std::max(t_1, t_2);
        } else {
            t = std::min(t_1, t_2);
        }

        intersect = Intersection((Shape *) &sphere, t);
        return true;
    }
}

// Function that (should have) determines if a pixel is in a shadow
// Method:
// dL = distance between ray origin and the light source or magnitude(lightPosition - shadowRayOrigin)
// T = distance between ray origin and the point that shadowRayNormal hits or magnitude(shadowRayOrigin + shadowRayNormal * it.tHit)
// if t < dL: pixel is in shadow of the light
bool inShadow(Tuple& shadowRayOrigin, Tuple& shadowRay, const std::vector<Node>& scene) {
    for (Node n : scene) {
        Intersection it;

        if (n.shapeType == 1) {
            Plane* obj = (Plane *) n.shapePtr;

            if (rayHitsPlane(shadowRayOrigin, shadowRay, *obj, it)) {
                if (it.tHit < 0) {
                    std::cout << "Plane t:" << it.tHit << std::endl;
                }
                return true;
            }
        } else {
            Sphere* obj = (Sphere *) n.shapePtr;

            if (rayHitsSphere(shadowRayOrigin, shadowRay, *obj, it)) {
                if (it.tHit < 0) {
                    std::cout << "Sphere t:" << it.tHit << std::endl;
                }
                return true;
            }
        }
    }

    return false;
}



void RayTrace(const std::vector<Node>& scene, const std::vector<LightSrc>& lSource) {
    double front_clip = 6.0;
    double w = 6.0;
    double h = 6.0;

    // Bottom-left point of image output window
    Tuple b_left = Tuple(-w/2, -h/2, front_clip, 1);

    // Point/Vector constants that will be reused in the main loop
    Tuple X = Tuple(1.0, 0, 0);
    Tuple Y = Tuple(0, 1.0, 0);
    Tuple cameraPoint = Tuple(0,0,0, 1);

    int image_pixel_size = 300;
    PPM ray_trace_image = easyppm_create(image_pixel_size, image_pixel_size, IMAGETYPE_PPM);
    easyppm_clear(&ray_trace_image, easyppm_rgb(255, 255, 255));

    for (double t = 0.0; t < h; t += (h / image_pixel_size)) {
        for (double s = 0.0; s < w; s += (w / image_pixel_size)) {
            // Calculate point P
            Tuple P = b_left + s * X + t * Y;

            // Create normalized ray from origin to P
            Tuple ray = P - cameraPoint;
            ray.normalize();

            // Declare intersection vector
            std::vector<Intersection> hits;

            // Iterate over every object in scene
            for (Node n : scene) {
                Intersection it;

                // Cast node to correct shape
                // shapeType == 1 -> plane
                // shapeType == 2 -> sphere
                if (n.shapeType == 1) {
                    Plane* obj = (Plane *) n.shapePtr;

                    if (rayHitsPlane(cameraPoint, ray, *obj, it)) {
                        hits.push_back(it);
                    }
                } else {
                    //std::cout << "Foo" << std::endl;
                    Sphere* obj = (Sphere *) n.shapePtr;

                    if (rayHitsSphere(cameraPoint, ray, *obj, it)) {
                        hits.push_back(it);
                    }
                }
            }

//            for (Node n : scene) {
//                std::cout << n.shapeType << std::endl;
//            }
//            std::cout << "-------" << std::endl;


            // Find the closest intersection, unpack shape and find point of intersection
            sort(hits.begin(), hits.end(), sortByHT);
            Intersection minHit = hits[0];
            Shape curObj = *minHit.obj;
            Tuple intersectPoint = cameraPoint + ray * minHit.tHit;

            // Calculate lighting from all light sources in the scene
            Tuple lightNormal;
            if (curObj.type == 1) {
                lightNormal = Tuple(0, 0, -1);
            } else {
                lightNormal = intersectPoint - curObj.position;
                lightNormal.normalize();
            }

//          This was my attempt at making shadows work, couldn't track down exactly why this never worked.
//          When used the entire image is shadowed i.e. ambient lighting is the only thing accounted for
//          Wanted to leave this in so you could see my attempt and maybe explain what was wrong in my feedback.
            Rgb totalLight = Rgb();
            for (LightSrc src : lSource) {
                // Handle shadows
                // Compute new point at an offset from intersectionPoint in the direction of the light source to avoid self-intersection
                Tuple shadowRay = src.position - intersectPoint;
                Tuple shadowRayOrigin = intersectPoint + shadowRay;
                bool hasShadow = inShadow(shadowRayOrigin, shadowRay, scene);

                // Compute ambient light regardless, only compute specular and diffuse lighting if point is not in shadows
                totalLight = totalLight + lightAmbient(curObj.rAmb, src.iAmb);
                if (not hasShadow) {
                    totalLight = totalLight +
                                 lightDiffuse(curObj.rDiff, curObj.position, lightNormal, src.iDiff, src.position);
                    totalLight = totalLight +
                                 lightSpecular(curObj.rSpec, curObj.position, lightNormal, src.iSpec, src.position,
                                               cameraPoint, src.specExp);
//                    std::cout << "Working" << std::endl;
                }
            }

//            Rgb totalLight = Rgb();
//            for (LightSrc src : lSource) {
//                // Compute ambient light regardless, only compute specular and diffuse lighting if point is not in shadows
//                totalLight = totalLight + lightAmbient(curObj.rAmb, src.iAmb);
//                totalLight = totalLight + lightDiffuse(curObj.rDiff, curObj.position, lightNormal, src.iDiff, src.position);
//                totalLight = totalLight + lightSpecular(curObj.rSpec, curObj.position, lightNormal, src.iSpec, src.position,
//                                               cameraPoint, src.specExp);
//            }


            int x = (int)((s * image_pixel_size / w) * (1.0 + std::numeric_limits<double>::epsilon()));
            int y = (int)((t * image_pixel_size / h) * (1.0 + std::numeric_limits<double>::epsilon()));
            easyppm_set(&ray_trace_image, x, y, easyppm_rgb(255 * totalLight.getR(),
                                                                    255 * totalLight.getG(),
                                                                    255 * totalLight.getB()));
        }
    }

//    for (Node n : scene) {
//        std::cout << n.shapeType << std::endl;
//    }

    // Write PPM image to file
    easyppm_write(&ray_trace_image, "image.ppm");
    easyppm_destroy(&ray_trace_image);
}


int main() {
    Plane p1 = Plane(Tuple(0, -4, 20, 1), Tuple(-1, 0, -1), Rgb(1, 0, 0),
                     Rgb(1, 0, 0), Rgb(1, 0, 0));
    Node o1 = {&p1, p1.type};

    Plane p2 = Plane(Tuple(0, -4, 20, 1), Tuple(1, 0, -1), Rgb(0, 0, 1),
                     Rgb(0, 0, 1), Rgb(0, 0, 1));
    Node o2 = {&p2, p2.type};

    Sphere s1 = Sphere(Tuple(1.5, 1.5, 6, 1), 1, Rgb(0, 255, 0),
                       Rgb(0, 1, 0), Rgb(0, 1, 0));
    Node o3 = {&s1, s1.type};

    Sphere s2 = Sphere(Tuple(0, 1, 10, 1), 2, Rgb(0.8,0.6,0.8),
                       Rgb(0.8,0.6,0.8), Rgb(0.8,0.6,0.8));
    Node o4 = {&s2, s2.type};

    Sphere s3 = Sphere(Tuple(-1.5, -1.5, 6, 1), 1, Rgb(0.5,0,0.5),
                       Rgb(0.5,0,0.5), Rgb(0.5,0,0.5));
    Node o5 = {&s3, s3.type};

    std::vector<Node> scene {o1, o2, o3, o4, o5};

    LightSrc l1 = {Tuple(3, 3, 6, 1),
                   Rgb(0.3, 0.3, 0.3),
                   Rgb(0.5, 0.5, 0.5),
                   Rgb(0.8, 0.8, 0.8), 2};
    LightSrc l2 = {Tuple(-3, 3, 6, 1),
                   Rgb(0.3, 0.3, 0.3),
                   Rgb(0.2, 0.2, 0.2),
                   Rgb(0.4, 0.4, 0.4), 1};

    std::vector<LightSrc> lights {l1, l2};

    RayTrace(scene, lights);

    return 0;
}
