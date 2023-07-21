#pragma  once
#include "../Common/math.hpp"
#include <optional>
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
    std::optional<vec3> tangent= std::nullopt;
    vec2 uv;
    Float epsilon = Constant::EPSILON;
    Spectrum Le(const vec3 & wo) const;
    Ray spawnRay(const vec3 & target){
        Float d = distance(p,target);
        vec3 dir = (target-p)/d;
        return Ray(p,dir,epsilon,d-epsilon);
    }
};
