#include "InfiniteSphereCap.h"
#include "scene.hpp"
static inline vec3 uniformSphericalCap(const vec2 &xi, float cosThetaMax)
{
    float phi = xi.x*Constant::TWO_PI;
    float z = xi.y*(1.0f - cosThetaMax) + cosThetaMax;
    float r = std::sqrt(std::max(1.0f - z*z, 0.0f));
    return vec3(
            std::cos(phi)*r,
            std::sin(phi)*r,
            z
    );
}

static inline float uniformSphericalCapPdf(float cosThetaMax)
{
    return  Constant::INV_TWO_PI/(1.0f - cosThetaMax);
}

Spectrum InfinteSphereCap::sampleLi(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf,
                                    VisibilityTester * vis) const {
   *wi = _capFrame.toWorld(uniformSphericalCap(u,_cosCapAngle));
   *pdf = uniformSphericalCapPdf(_cosCapAngle);

   Intersection pShape;
   pShape.p = ref.p + _worldRadius * 2 * *wi;
    * vis = VisibilityTester(ref, pShape);

   return _emission;
}

LightSampleResult InfinteSphereCap::sampleDirect(const vec2 & positionSample, const vec2 & u2) {
    _NOT_IMPLEMENT_ERROR;
}

Spectrum InfinteSphereCap::Le(const Ray & ray) const {

    if ( dot(ray.d,_capDir) < _cosCapAngle)
        return Spectrum(0);
    return _emission ;
}

Float InfinteSphereCap::PdfLi(const Intersection & pShape, const vec3 & ref) const {
    vec3 rayDir = -pShape.w;
    if(dot(rayDir,_capDir)<_cosCapAngle)
        return 0;
    return uniformSphericalCapPdf(_cosCapAngle);
}

void InfinteSphereCap::logDebugInfo( ) const {

}

Spectrum InfinteSphereCap::Power( ) {
    return 2 * Constant::PI * _emission;
}

void InfinteSphereCap::Preprocess(const Scene & scene) {
    scene.getWorldBound().BoundingSphere(& _worldCenter, & _worldRadius);
}

InfinteSphereCap::InfinteSphereCap(const Json & json) : Light(int(LightFlags::Infinite)) {

    _capAngle = Angle::degToRad(getOptional(json,"cap_angle",5));
    _cosCapAngle = cos(_capAngle);

    mat4 transform = getOptional(json,"transform",getIndentifyTransform());
    _capDir = transformVector(transform,vec3(0,1,0));
    _capFrame = Frame(_capDir);

    _emission = getOptional(json,"emission",Spectrum(1));
    if(json.contains("power")){
        Spectrum  power = json["power"];
        _emission = power * Constant::INV_TWO_PI/(1.0f - _cosCapAngle);
    }
}