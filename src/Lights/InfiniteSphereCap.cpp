#include "InfiniteSphereCap.h"
#include "scene.hpp"
#include "SampleRecords/PositionAndDirectionSample.h"
#include "Sampler/Warp.hpp"

static inline vec3 uniformSphericalCap(const vec2 & xi, float cosThetaMax) {
    float phi = xi.x * Constant::TWO_PI;
    float z = xi.y * ( 1.0f - cosThetaMax ) + cosThetaMax;
    float r = std::sqrt(std::max(1.0f - z * z, 0.0f));
    return vec3(
            std::cos(phi) * r,
            std::sin(phi) * r,
            z
    );
}

static inline float uniformSphericalCapPdf(float cosThetaMax) {
    return Constant::INV_TWO_PI / ( 1.0f - cosThetaMax );
}

Spectrum InfinteSphereCap::sampleLi(const vec3 & ref, const vec2 & u, vec3 * wi, Float * pdf,
                                    Float * distance) const {
    * wi = _capFrame.toWorld(uniformSphericalCap(u, _cosCapAngle));
    * pdf = uniformSphericalCapPdf(_cosCapAngle);
    * distance = 2 * _worldRadius;
    return _emission;
}

PositionAndDirectionSample InfinteSphereCap::sampleDirect(const vec2 & positionSample, const vec2 & dirSample) const{
    PositionAndDirectionSample sample;
    auto dir = _capFrame.toWorld(uniformSphericalCap(dirSample,_cosCapAngle));
    sample.dirPdf = uniformSphericalCapPdf(_cosCapAngle);
    sample.n = dir;
    vec3 v1, v2;
    coordinateSystem(- dir, v1, v2);
    vec2 cd = Warp::ConcentricSampleDisk(positionSample);
    vec3 pDisk = _worldCenter + ( cd.x * v1 + cd.y * v2 ) * _worldRadius;
    sample.posPdf  =1 / (_worldRadius * _worldRadius * Constant::PI );
    sample.weight =  _emission;
    _NOT_IMPLEMENT_ERROR;
}

Spectrum InfinteSphereCap::Le(const Ray & ray) const {

    if ( dot(ray.d, _capDir) < _cosCapAngle )
        return Spectrum(0);
    return _emission;
}

Float InfinteSphereCap::PdfLi(const Intersection & pShape, const vec3 & ref) const {
    vec3 rayDir =  pShape.w;
    if ( dot(rayDir, _capDir) < _cosCapAngle )
        return 0;
    return uniformSphericalCapPdf(_cosCapAngle);
}

void InfinteSphereCap::logDebugInfo( ) const {

}

Spectrum InfinteSphereCap::Power( ) {
    return 2 * Constant::PI * _emission;
}


InfinteSphereCap::InfinteSphereCap(const Json & json)  {

    _capAngle = Angle::degToRad(getOptional(json, "cap_angle", 5));
    _cosCapAngle = abs(cos(_capAngle));
    _cosCapAngle = std::min(_cosCapAngle,0.99f);
    mat4 transform = getOptional(json, "transform", getIndentifyTransform());
    _capDir = transformVector(transform, vec3(0, 1, 0));
    _capFrame = Frame(_capDir);

    _emission = getOptional(json, "emission", Spectrum(1));
    if ( json.contains("power") ) {
        Spectrum power = json["power"];
        _emission = power * Constant::INV_TWO_PI / ( 1.0f - _cosCapAngle );
    }
}
