// CS360 Project 5 - Gabe Kelly
#include <vector>
#include <limits>
#include <algorithm>

#include "Lighting.h"
#include "Intersection.h"
#include "Plane.h"
#include "Sphere.h"
#include "easyppm.h"

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
                break;
            }
        } else {
            auto obj = (Sphere *) n.shapePtr;

            if (rayHitsSphere(intersectPoint, intersectToLight, *obj, it)) {
                t = it.tHit;
                break;
            }
        }
    }

    return t < dL;
}

Rgb Trace(Ray& ray, const std::vector<Node>& scene, const std::vector<LightSrc>& lSource) {
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


    Rgb totalLight = Rgb();
    for (LightSrc src : lSource) {
        // Compute ambient light regardless, only compute specular and diffuse lighting if point is not in shadows
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

    return totalLight;
}

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
            Rgb totalLight = Trace(ray, scene, lSource);

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

    Render(scene, lights);

    return 0;
}
