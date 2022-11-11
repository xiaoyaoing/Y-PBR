#include "Distant.hpp"
#include "scene.hpp"
Float DistantLight::PdfLi(const Intersection & pShape, const vec3 & ref) const {
    return 0;
}

Spectrum
DistantLight::sampleLi(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf, VisibilityTester * vis) const {
    *wi = wLight;
    *pdf = 1;
    vec3 pOutside = ref.p + wLight * (2 * worldRadius);
    Intersection pShape;
    pShape.p = pOutside;
    *vis = VisibilityTester(ref,pShape);
    return L;
}

LightSampleResult DistantLight::sampleDirect(const vec2 & positionSample, const vec2 & dirSample) {
    _NOT_IMPLEMENT_ERROR;
}

Spectrum DistantLight::Le(const Ray & ray) const {
    return Spectrum(0);
}

Spectrum DistantLight::Power( ) {
    return L * Constant::PI * worldRadius * worldRadius;
}

void DistantLight::Preprocess(const Scene & scene) {
    scene.getWorldBound().BoundingSphere(&worldCenter, &worldRadius);
}

DistantLight::DistantLight(const Json & json) : Light(int(LightFlags::DeltaDirection)) {
    L = getOptional(json,"radiance",Spectrum(1));
    vec3 from = getOptional(json,"from",vec3(0,0,0));
    vec3 to = getOptional(json,"to",vec3(0,0,1));
    wLight = normalize((to-from));
}

bool DistantLight::isDeltaLight( ) const {
    return true;
}
