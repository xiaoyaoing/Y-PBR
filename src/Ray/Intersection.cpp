#include "Intersection.hpp"
#include "../Primitives/Primitive.hpp"

Spectrum Intersection::Le(const vec3& wo) const {
    const AreaLight* area = primitive->areaLight.get();
    return area ? area->directLighting(*this, wo) : Spectrum(0.f);
}