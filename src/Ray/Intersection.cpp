#include "Intersection.hpp"
#include "../Primitives/Primitive.hpp"

Spectrum Intersection::Le(const vec3 & w) const {
    const AreaLight *area = primitive->areaLight.get();
    return area ? area->L(*this, w) : Spectrum(0.f);
}

vec3 Intersection::toLocal(const vec3 & w) const {
    return normalize(shFrame.toLocal(w));
}

vec3 Intersection::toWorld(const vec3 & w) const  {
    return normalize(shFrame.toWorld(w));
}

void Intersection::setNormal(const vec3 & normal)  {
    this->n=normal;
    this->shFrame=Frame(n);
}
