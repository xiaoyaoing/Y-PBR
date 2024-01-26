#include "Distant.hpp"
#include "scene.hpp"
#include "SampleRecords/PositionAndDirectionSample.h"

Float DistantLight::PdfLi(const Intersection& pShape, const vec3& ref) const {
    return 0;
}

Spectrum
DistantLight::sampleLi(const vec3& ref, const vec2& u, vec3* wi, Float* pdf, Float* distance) const {
    *wi       = wLight;
    *pdf      = 1;
    *distance = 2 * _worldRadius;
    return L;
}

PositionAndDirectionSample DistantLight::sampleDirect(const vec2& positionSample, const vec2& dirSample) const {
    _NOT_IMPLEMENT_ERROR;
}

Spectrum DistantLight::Le(const Ray& ray) const {
    return Spectrum(0);
}

Spectrum DistantLight::Power() {
    return L * Constant::PI * _worldRadius * _worldRadius;
}

DistantLight::DistantLight(const Json& json) : Infinite(LightFlags::DeltaDirection) {
    L         = getOptional(json, "radiance", Spectrum(1));
    vec3 from = getOptional(json, "from", vec3(0, 0, 0));
    vec3 to   = getOptional(json, "to", vec3(0, 0, 1));
    wLight    = normalize((to - from));
}

bool DistantLight::isDeltaLight() const {
    return true;
}