// CS360 Project 5 - Gabe Kelly
#include <vector>
#include <limits>
#include <algorithm>
#include <string>

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
    Shape *shapePtr;
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

bool rayHitsPlane(const Tuple &rayOriginPoint, const Tuple &rayDirectionVector, const Plane &plane,
                  Intersection &intersect) {
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

bool rayHitsSphere(const Tuple &rayOriginPoint, const Tuple &rayDirectionVector, const Sphere &sphere,
                   Intersection &intersect) {
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

        // Two positive roots
        if (t_1 > 0 and t_2 > 0) {
            t = min(t_1, t_2);
        }
            // Only 1 positive root
        else if ((t_1 < 0) != (t_2 < 0)) {
            t = max(t_1, t_2);
        }
            // No positive roots
        else {
            return false;
        }

        intersect = Intersection((Shape *) &sphere, t);
        return true;
    }
}

bool inShadow(Tuple intersectPoint, Tuple &lightPoint, const std::vector<Node> &scene) {
    intersectPoint = intersectPoint + (0.001 * intersectPoint - lightPoint);
    Tuple intersectToLight = lightPoint - intersectPoint;
    double dL = intersectToLight.magnitude();
    intersectToLight.normalize();

    Intersection it;
    double t;
    for (Node n: scene) {
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

Rgb Trace(Ray &ray, const std::vector<Node> &scene, const std::vector<LightSrc> &lSource, int recursion_depth) {
    // Declare intersection vector
    std::vector<Intersection> hits;

    // Iterate over every object in scene
    for (Node n: scene) {
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

    // If nothing is hit, return background color
    if (hits.empty()) {
        return {0.856, 0.952, 0.992};
    }

    // Find the closest intersection, unpack shape and find point of intersection
    sort(hits.begin(), hits.end(), sortByHT);

    // Get hit object
    Intersection minHit = hits[0];
    Shape curObj = *minHit.obj;
    Tuple intersectPoint = ray.point + ray.direction * minHit.tHit;

    // Calculate lighting from all light sources in the scene
    Tuple hitNormal;
    if (curObj.type == 1) {
        hitNormal = Tuple(0, 0, -1);

    } else {
        hitNormal = intersectPoint - curObj.position;
        hitNormal.normalize();

//        if (ray.direction.dot(hitNormal) < 0) {
//            hitNormal.z *= -1;
//        }
        if ((ray.direction.dot(hitNormal) < 0) and (curObj.reflectiveness >= 1)) {
            hitNormal.z *= -1;
        }
    }

    Rgb totalLight = Rgb();
    for (LightSrc src: lSource) {
        totalLight = totalLight + lightAmbient(curObj.rAmb, src.iAmb);
        bool hasShadow = inShadow(intersectPoint, src.position, scene);

        if (not hasShadow and curObj.reflectiveness >= 1) {
            totalLight = totalLight +
                         lightDiffuse(curObj.rDiff, curObj.position, hitNormal, src.iDiff, src.position);
            totalLight = totalLight +
                         lightSpecular(curObj.rSpec, curObj.position, hitNormal, src.iSpec, src.position,
                                       ray.point, src.specExp);
        }
    }

    if ((curObj.reflectiveness < 1) and (recursion_depth < MAX_RECURSION_DEPTH)) {
        // Compute reflected ray
        Tuple reflectDirection = -2 * ray.direction.dot(hitNormal) * hitNormal - ray.direction;
        reflectDirection.normalize();

        Tuple reflectOrigin = intersectPoint + (hitNormal * 0.001);
        Ray reflectionRay = {reflectOrigin, reflectDirection};
        totalLight = totalLight + (curObj.reflectiveness * Trace(reflectionRay, scene, lSource, recursion_depth + 1));
    }

    return totalLight;
}

void Render(const std::vector<Node> &scene, const std::vector<LightSrc> &lSource, std::string &outputFileName) {
    double front_clip = 6.0;
    double w = 6.0;
    double h = 6.0;

    // Bottom-left point of image output window
    Tuple b_left = Tuple(-w / 2, -h / 2, front_clip, 1);

    // Point/Vector constants that will be reused in the main loop
    Tuple X = Tuple(1.0, 0, 0);
    Tuple Y = Tuple(0, 1.0, 0);
    Tuple cameraPoint = Tuple(0, 0, 0, 1);

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
            int x = (int) ((s * image_pixel_size / w) * (1.0 + std::numeric_limits<double>::epsilon()));
            int y = (int) ((t * image_pixel_size / h) * (1.0 + std::numeric_limits<double>::epsilon()));
            easyppm_set(&ray_trace_image, x, y, easyppm_rgb(255 * totalLight.getR(),
                                                            255 * totalLight.getG(),
                                                            255 * totalLight.getB()));
        }
    }

    // Write PPM image to file
    easyppm_write(&ray_trace_image, outputFileName.c_str());
    easyppm_destroy(&ray_trace_image);
}


void GenerateWaveAnimation() {
    LightSrc l1 = {Tuple(3, 3, 6, 1),
                   Rgb(0.3, 0.3, 0.3),
                   Rgb(0.5, 0.5, 0.5),
                   Rgb(0.8, 0.8, 0.8), 2};
    LightSrc l2 = {Tuple(-3, 3, 6, 1),
                   Rgb(0.3, 0.3, 0.3),
                   Rgb(0.2, 0.2, 0.2),
                   Rgb(0.4, 0.4, 0.4), 1};
    const std::vector<LightSrc> lights{l1, l2};

    Sphere s1 = Sphere(Tuple(0, 1, 30, 1), 3, 0.4, Rgb(0.8, 0.6, 0.8),
                       Rgb(0.8, 0.6, 0.8), Rgb(0.8, 0.6, 0.8));
    Node o1 = {&s1, s1.type};

    Sphere s2 = Sphere(Tuple(-2.5, 5, 20, 1), 1, 1, Rgb(0.5, 0, 0.5),
                       Rgb(0.5, 0, 0.5), Rgb(0.5, 0, 0.5));
    Node o2 = {&s2, s2.type};

    Sphere s3 = Sphere(Tuple(2.5, -5, 20, 1), 1, 1, Rgb(0, 0.5, 0.5),
                       Rgb(0, 0.5, 0.5), Rgb(0, 0.5, 0.5));
    Node o3 = {&s3, s3.type};

    Sphere s4 = Sphere(Tuple(5.0, 5.0, 25, 1), 2, 1, Rgb(0.2, 0.6, 0.8),
                       Rgb(0.2, 0.6, 0.8), Rgb(0.2, 0.6, 0.8));
    Node o4 = {&s4, s4.type};

    std::vector<Node> spheres{o1, o2, o3, o4};

    int frames = 60;
    double t_max = 6.3;

    // Loop to create each frame
    for (int iter = 0; iter <= frames; iter++) {
        // Calculate delta-value
        double delta = iter * (t_max / frames);
        for (int i = 1; i <= spheres.size() - 1; i++) {
            Tuple displacement = Tuple(4 * cos(delta + i),
                                       16 * sin(delta + i), 3);
            displacement.normalize();
            spheres[i].shapePtr->position = spheres[i].shapePtr->position + displacement;

        }

        std::string fileName = "animation_frames/frame_X.ppm";
        size_t replace_index = "animation_frames/frame_X.ppm"sv.find('X');
        fileName.replace(replace_index, 1, std::to_string(iter));
        Render(spheres, lights, fileName);
        std::cout << "Rendering frame..." << std::endl;
    }

    std::cout << "Done!" << std::endl;

}

int main() {
    Sphere s1 = Sphere(Tuple(1.5, 1.5, 6, 1), 1, 1, Rgb(0, 1, 0),
                       Rgb(0, 1, 0), Rgb(0, 1, 0));
    Node o3 = {&s1, s1.type};

    Sphere s2 = Sphere(Tuple(0, 1, 10, 1), 2, 0.4, Rgb(0.8, 0.6, 0.8),
                       Rgb(0.8, 0.6, 0.8), Rgb(0.8, 0.6, 0.8));
    Node o4 = {&s2, s2.type};

    Sphere s3 = Sphere(Tuple(-1.5, -1.5, 6, 1), 1, 1, Rgb(0.5, 0, 0.5),
                       Rgb(0.5, 0, 0.5), Rgb(0.5, 0, 0.5));
    Node o5 = {&s3, s3.type};

    Sphere s4 = Sphere(Tuple(-1.5, 1.5, 6, 1), 1, 1, Rgb(0, 0.5, 0.5),
                       Rgb(0, 0.5, 0.5), Rgb(0, 0.5, 0.5));
    Node o6 = {&s4, s4.type};

    Sphere s5 = Sphere(Tuple(2.0, -1.5, 14, 1), 2, 1, Rgb(0.2, 0.6, 0.8),
                       Rgb(0.2, 0.6, 0.8), Rgb(0.2, 0.6, 0.8));
    Node o7 = {&s5, s5.type};

    std::vector<Node> scene{o3, o4, o5, o6, o7};

    LightSrc l1 = {Tuple(3, 3, 6, 1),
                   Rgb(0.3, 0.3, 0.3),
                   Rgb(0.5, 0.5, 0.5),
                   Rgb(0.8, 0.8, 0.8), 2};
    LightSrc l2 = {Tuple(-3, 3, 6, 1),
                   Rgb(0.3, 0.3, 0.3),
                   Rgb(0.2, 0.2, 0.2),
                   Rgb(0.4, 0.4, 0.4), 1};

    std::vector<LightSrc> lights{l1, l2};

    std::string fileName = "image.ppm";
    Render(scene, lights, fileName);

//    GenerateWaveAnimation();

    return 0;
}