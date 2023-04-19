#pragma  once
#include "Ray/Ray.hpp"
#include "Colors/Spectrum.hpp"
struct PositionAndDirectionSample{
    vec3 n;
    Float posPdf;
    Float dirPdf;
    Spectrum weight;
    Ray ray;
    Spectrum computeWeight(){
        return radiance * absDot(ray.d,n) / (posPdf * dirPdf);
    }
};

