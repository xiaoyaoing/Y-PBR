#pragma  once
#include "../Common/math.hpp"
#include "../Bsdfs/Reflection.hpp"
#include "../Common/Frame.hpp"
#include "../Colors/Spectrum.hpp"

class Primitive ;

struct Intersection {
    vec3 p;                    //position
    vec3 wo;                   //-ray.d(world space)
    Bsdf * bsdf;
    const Primitive * primitive;// Primitive ptr

    void setNormal(const vec3 & n); //set normalAndFrame
    vec3 getNormal() const  {return n; }

    vec3 toLocal(const vec3 & w) const ;

    vec3 toWorld(const vec3 & w) const ;

    Spectrum Le(const vec3 & w) const;

private:
    vec3 n;                    //normal
    Frame shFrame;             //shading Frame

};