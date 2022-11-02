#pragma  once
#include "../Common/math.hpp"
//#include "../Bsdfs/Reflection.hpp"

#include "../Colors/Spectrum.hpp"

class Primitive ;
class BSDF;


struct Intersection {
    vec3 p;
    vec3 w;  //dir of the indicent ray
    BSDF * bsdf;
    const Primitive * primitive;
    vec3 Ns ;                    //Geometric Normal
    vec3 Ng ;                    //Shading Normal
    vec2 uv;

    Spectrum Le(const vec3 & w) const;
    Intersection(){}
    Intersection(const Intersection & its):p(its.p),w(its.w),bsdf(its.bsdf),Ns(its.Ns),Ng(its.Ng),uv(its.uv),primitive(its.primitive){

    }
};

struct SurfaceIntersection : Intersection{

};
