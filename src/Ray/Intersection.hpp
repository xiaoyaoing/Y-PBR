#pragma  once
#include "../Common/math.hpp"
//#include "../Bsdfs/Reflection.hpp"

#include "../Colors/Spectrum.hpp"
#include "Ray.hpp"

class Primitive ;
class BSDF;
class BSSRDF;

struct Intersection {
    vec3 p;
    vec3 w;  //dir of the indicent ray
    BSDF * bsdf = nullptr;
    BSSRDF * bssrdf = nullptr;
    const Primitive * primitive;
    vec3 Ns,Ng ;
    vec3 * tangent = nullptr;
    vec2 uv;
    Float epsilon = Constant::EPSILON;
    Spectrum Le(const vec3 & wo) const;
    Intersection(){}
    Intersection(const Intersection & its):p(its.p),w(its.w),bsdf(its.bsdf),bssrdf(its.bssrdf)
                                        ,Ns(its.Ns),Ng(its.Ng),uv(its.uv),primitive(its.primitive),tangent(its.tangent){

    }
    Ray spawnRay(const vec3 & target){
        Float d = distance(p,target);
        vec3 dir = (target-p)/d;
        return Ray(p,dir,epsilon,d-epsilon);
    }
};

struct SurfaceIntersection : Intersection{

};
