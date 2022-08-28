#include "Light.hpp"
#include "../Primitives/Primitive.hpp"
#include "../scene.hpp"


Spectrum
AreaLight::Sample_Li(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf, VisibilityTester * vis) const {

    Intersection pShape = primitive->Sample(ref,u,pdf);

    if(*pdf==0){
        return Spectrum();
    }

    *wi = normalize(ref.p-pShape.p);
    *vis = VisibilityTester(ref, pShape);

    auto dotN = dot(*wi,pShape.getNormal());
    if(ref.bsdf->name=="shortBox")
        int k=1;
    return L(pShape,*wi);

}

Spectrum AreaLight::L(const Intersection & intr, const vec3 & w) const {
    return (twoSide || dot(intr.getNormal(), w) > 0) ? albedo : Spectrum(0.f);
}



AreaLight::AreaLight(const std::shared_ptr< Primitive > & primitive,
                     const Spectrum & albedo,
                     bool twoSide)
                        :primitive(primitive),albedo(albedo),twoSide(twoSide)
{

}

bool VisibilityTester::Unoccluded(const Scene & scene) const {
    vec3 dir =(p1.p-p0.p);
    Float distance = length(dir);
    dir= normalize(dir);
    Ray ray(p0.p+Constant::EPSILON * dir ,dir,Constant::EPSILON,distance-10*Constant::EPSILON);
    return !scene.intersectP(ray);

}
