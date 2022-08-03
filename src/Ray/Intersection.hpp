#pragma  once
#include "../Common/math.hpp"
#include "../Bsdfs/Bsdf.hpp"
#include "../Common/Frame.hpp"
#include "../Colors/Spectrum.hpp"
struct Intersection {
    vec3 p;
    vec3 n;
    vec3 wo;
    Bsdf * bsdf ;
    Frame shFrame;

    Spectrum Le(const vec3 & w);


};