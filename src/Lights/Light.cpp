#include "Light.hpp"
#include "../Primitives/Primitive.hpp"

Spectrum
AreaLight::Sample_Li(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf, VisibilityTester * vis) const {

    Intersection pShape = primitive->Sample(ref,u,pdf);

    *wi = normalize(ref.p-pShape.p);
    *vis = VisibilityTester(ref, pShape);

    return L(pShape,*wi);

}

Spectrum AreaLight::L(const Intersection & intr, const vec3 & w) const {
    return (twoSide || dot(intr.n, w) > 0) ? albedo : Spectrum(0.f);
}



AreaLight::AreaLight(const std::shared_ptr< Primitive > & primitive,
                     const Spectrum & albedo,
                     bool twoSide)
                        :primitive(primitive),albedo(albedo),twoSide(twoSide)
{

}
