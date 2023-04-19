#include "InfiniteSphere.hpp"
#include "scene.hpp"
#include "iostream"
#include "Sampler/Warp.hpp"
#include <spdlog/spdlog.h>

#include "Common/Transform.hpp"
#include "SampleRecords/PositionAndDirectionSample.h"

static vec3 rgb(0);
static int count = 0;

Spectrum InfinteSphere::Le(const Ray & ray) const {
    vec2 uv = directionToUV(ray.d);
    //return Spectrum(0,uv.y,0);
    //return (ray.d + 1.0f)/2.f;
    Spectrum L = _emission->eval(uv);
    if ( hasNan(L) ) {
        uv = directionToUV(ray.d);
        L = _emission->eval(uv);
    }
    return L;
}

Spectrum InfinteSphere::sampleLi(const vec3 & ref, const vec2 & u, vec3 * wi, Float * pdf,
                                 Float * distance) const {

//    auto t = getTransFormMfdiratrix(vec3(0),vec3(1),vec3(270,270,270)) * this->_toWorld;
//    auto s = Mat4ToStr(t);
    Float mapPdf;
    vec2 uv = _emission->sample(MAP_SPHERICAL, u, & mapPdf);
    if ( ! mapPdf )
        return Spectrum(0);
    float sinTheta;
    * wi = uvToDirection(uv, sinTheta);
    if ( sinTheta == 0 )
        * pdf = 0;
    else
        * pdf = mapPdf / ( 2 * Constant::PI * Constant::PI * sinTheta );
    * distance = 2 * _worldRadius;
    //return Spectrum(2);
    return _emission->eval(uv);
}

//Spectrum InfinteSphere::directLighting(const Intersection & intr) const {
//    return _emission->eval(directionToUV(intr.w));
//}

Spectrum InfinteSphere::Power( ) {
    return 4 * Constant::PI * _worldRadius * _worldRadius * _emission->average();
}

void InfinteSphere::Preprocess(const Scene & scene) {
    Infinite::Preprocess(scene);
    this->_emission->makeSamplable(MAP_SPHERICAL);
}

vec2 InfinteSphere::directionToUV(const vec3 & wi) const {
    vec3 wLocal = transformVector(_toLocal, wi);
    auto t = length(wi);
    auto t1 = length(wLocal);
    if ( abs(wLocal.y) > 1 ) {

    }
    auto uv = vec2(std::atan2(wLocal.z, wLocal.x) * Constant::INV_TWO_PI + 0.5f,
                   std::acos(- clamp(wLocal.y, -1, 1)) * Constant::INV_PI);
    return uv;
}

vec2 InfinteSphere::directionToUV(const vec3 & wi, float & sinTheta) const {
    vec3 wLocal = transformVector(_toLocal, wi);
    sinTheta = sqrt(1 - wLocal.y * wLocal.y);
    return vec2(std::atan2(wLocal.z, wLocal.x) * Constant::INV_TWO_PI + 0.5f,
                std::acos(- clamp(wLocal.y, -1, 1)) * Constant::INV_PI);
}

vec3 InfinteSphere::uvToDirection(vec2 uv, float & sinTheta) const {
    float phi = ( uv.x - 0.5f ) * 2 * Constant::PI;
    float theta = uv.y * Constant::PI;
    sinTheta = std::sin(theta);
    vec3 wLocal = vec3(
            std::cos(phi) * sinTheta,
            - std::cos(theta),
            std::sin(phi) * sinTheta);
    return transformVector(_toWorld, wLocal);
}

Float InfinteSphere::PdfLi(const Intersection & pShape, const vec3 /*ref*/&) const {
    Float sinTheta;
    vec2 uv = directionToUV(pShape.w, sinTheta);
    if ( sinTheta == 0 ) return 0;
    return Constant::INV_PI * Constant::INV_TWO_PI * _emission->pdf(MAP_SPHERICAL, uv) / sinTheta;
}

void InfinteSphere::logDebugInfo( ) const {
    spdlog::info("infinite" + toColorStr(rgb / (float) count));
}

PositionAndDirectionSample InfinteSphere::sampleDirect(const vec2 & positionSample, const vec2 & dirSample) const{
    PositionAndDirectionSample result;
    Float mapPdf;
    vec2 uv = _emission->sample(MAP_SPHERICAL, dirSample, & mapPdf);
    if ( ! mapPdf )
        return result;
    float sinTheta;
    vec3 dir = - uvToDirection(uv, sinTheta);
    result.n = dir;
//    *pdf = mapPdf / ( 2 * Constant::PI * Constant::PI * sinTheta );
    result.dirPdf = _emission->pdf(MAP_SPHERICAL, uv) / (2 * Constant::PI * Constant::PI * sinTheta );
    result.posPdf = 1 / (_worldRadius * _worldRadius * Constant::PI );
    result.weight = _emission->eval(uv) / (result.dirPdf * result.posPdf );
    vec3 v1, v2;
    coordinateSystem(- dir, v1, v2);
    vec2 cd = Warp::ConcentricSampleDisk(positionSample);
    vec3 pDisk = _worldCenter + ( cd.x * v1 + cd.y * v2 ) * _worldRadius;
    result.ray = Ray(pDisk + _worldRadius * - dir, dir);
    return result;
}

InfinteSphere::InfinteSphere(const Json & json) {
    _emission = std::make_shared < BitMapTexture < Spectrum > >
            (json.at("emission").get < std::string >());
    _emission->LoadResources();
    _toWorld = getOptional(json, "transform", getIndentifyTransform());
    _toLocal = glm::transpose(_toWorld);
}

