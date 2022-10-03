#include "Infinte.hpp"
#include "scene.hpp"
#include "iostream"
#include <spdlog/spdlog.h>

static vec3 rgb(0);
static int count=0;

Spectrum InfinteSphere::environmentLighting(const Ray & ray) const {
    return _emission->Evaluate(directionToUV(ray.d));
}

Spectrum InfinteSphere::Sample_Li(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf,
                                  VisibilityTester * vis) const {

    vec2 uv = _emission->sample(MAP_SPHERICAL, u);
    float sinTheta;
    *wi = uvToDirection(uv, sinTheta);
    * pdf = _emission->pdf(MAP_SPHERICAL, uv) / ( 2 * Constant::PI * Constant::PI * sinTheta );

    Intersection pShape;
    pShape.p = ref.p + 2 * _worldRadius * ( * wi );
    * vis = VisibilityTester(ref, pShape);

    return _emission->Evaluate(uv);
}

Spectrum InfinteSphere::directLighting(const Intersection & intr) const {
    return _emission->Evaluate(directionToUV(intr.w));
}

Float InfinteSphere::Power( ) {
    return Light::Power();
}

void InfinteSphere::Preprocess(const Scene & scene) {
    scene.getWorldBound().BoundingSphere(& _worldCenter, & _worldRadius);
    this->_emission->LoadResources();
    this->_emission->makeSamplable(MAP_SPHERICAL);
}

vec2 InfinteSphere::directionToUV(const vec3 & wi) const {
    vec3  wLocal = wi;
    return vec2 (std::atan2(wLocal.z, wLocal.x)*Constant::INV_TWO_PI+ 0.5f, std::acos(-wLocal.y)*Constant::INV_PI);
}

vec2 InfinteSphere::directionToUV(const vec3 & wi, float & sinTheta) const {
    vec3  wLocal = wi;
    sinTheta = sqrt(1-wLocal.y*wLocal.y);
    return vec2 (std::atan2(wLocal.z, wLocal.x)*Constant::INV_TWO_PI+ 0.5f, std::acos(-wLocal.y)*Constant::INV_PI);
}

vec3 InfinteSphere::uvToDirection(vec2 uv, float & sinTheta) const {
    float phi = ( uv.x-0.5f ) * 2*Constant::PI;
    float theta = uv.y * Constant::PI;
    sinTheta = std::sin(theta);
    return vec3(
//            std::cos(phi) * sinTheta,
//            - std::cos(theta),
//            std::sin(phi) * sinTheta
            std::cos(phi)*sinTheta,
            -std::cos(theta),
            std::sin(phi)*sinTheta
    );
}

Float InfinteSphere::directPdf(const Intersection & pShape, const vec3 /*ref*/&) const {

    Float sinTheta;
    vec2 uv = directionToUV(pShape.w, sinTheta);
    return Constant::INV_PI*Constant::INV_TWO_PI*_emission->pdf(MAP_SPHERICAL, uv)/sinTheta;

//    vec2 uv = directionToUV(- pShape.w);
//    return _emission->pdf(MAP_SPHERICAL, uv);
}

void InfinteSphere::logDebugInfo( ) const {
    {_emission->debugLoginfo();}
    spdlog::info("infinite"+toColorStr(rgb/(float)count));
}

