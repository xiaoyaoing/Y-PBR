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
    vec3 Ns,Ng ;
    vec3 tangent;
    vec2 uv;
    Float epsilon = Constant::EPSILON;
    Spectrum Le(const vec3 & wo) const;
    Intersection(){}
    Intersection(const Intersection & its):p(its.p),w(its.w),bsdf(its.bsdf),Ns(its.Ns),Ng(its.Ng),uv(its.uv),primitive(its.primitive),
    tangent(its.tangent){

    }
};

struct SurfaceIntersection : Intersection{

};
