#include "Intersection.hpp"
#include "../Primitives/Primitive.hpp"

Spectrum Intersection::Le(const vec3 & w) const {
    const AreaLight * area = primitive->areaLight.get();
    return area ? area->directLighting(* this) : Spectrum(0.f);


}