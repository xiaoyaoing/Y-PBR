#pragma  once
#include "Ray/Ray.hpp"
#include "Colors/Spectrum.hpp"
struct PositionAndDirectionSample{
    vec3 n;
    Float posPdf = 0;
    Float dirPdf = 0;
    Spectrum weight = Spectrum(0);
    Ray ray;
    Spectrum computeWeight(){
        return weight * absDot(n,ray.d) /(dirPdf * posPdf);
    }
};

