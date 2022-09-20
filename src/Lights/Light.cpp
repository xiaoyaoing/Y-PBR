#include "Light.hpp"
#include "../Primitives/Primitive.hpp"
#include "../scene.hpp"


Spectrum
AreaLight::Sample_Li(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf, VisibilityTester * vis) const {

    Intersection pShape = primitive->Sample(ref,u,pdf);

    if(*pdf==0){
        return Spectrum();
    }

    *wi = normalize(pShape.p-ref.p);
    pShape.w = *wi;
    *vis = VisibilityTester(ref, pShape);

    return directLighting(pShape);

}



Spectrum AreaLight::directLighting(const Intersection & intr) const {
    return (twoSide || dot(intr.Ng, -intr.w) > 0) ? albedo : Spectrum(0.f);
}



AreaLight::AreaLight(const std::shared_ptr< Primitive > & primitive,
                     const Spectrum & albedo,
                     bool twoSide)
                        :Light((int)LightFlags::Area),
                        primitive(primitive),albedo(albedo),twoSide(twoSide)
{

}

Float AreaLight::directPdf(const Intersection & pShape, const vec3 & ref) const {
    return primitive->directPdf(pShape,ref);
}

bool VisibilityTester::Unoccluded(const Scene & scene) const {
    vec3 dir =(p1.p-p0.p);
    Float distance = length(dir);
    dir= normalize(dir);
    Ray ray(p0.p ,dir,Constant::EPSILON,distance-Constant::EPSILON);
    return !scene.intersectP(ray);
}
