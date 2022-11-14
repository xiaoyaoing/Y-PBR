//2022/7/17

#include "Sphere.hpp"
#include "../Sampler/Warp.hpp"
#include "Common/Transform.hpp"

std::optional < Intersection > Sphere::intersect(Ray & ray) const {
    Intersection intersection;
    vec3 p = ray.o - center;
    float B = dot(p, ray.d);
    float C = length2(p) - pow2(radius);
    float detSq = B * B - C;
    if ( detSq >= 0.0f ) {
        float det = std::sqrt(detSq);
        float t = - B - det;

        if ( t > ray.farT || t  < ray.nearT ) {
            t+=2 * det;
            if(t > ray.farT || t  < ray.nearT)
            return std::nullopt;
        }
        ray.farT=t;
        intersection.p = ray(t);
        intersection.Ns = intersection.Ng = normal(intersection.p);
        intersection.primitive = this;
        intersection.bsdf = bsdf.get();
        intersection.w = ray.d;
        return {intersection};
    }
    return std::nullopt;
}

bool Sphere::occluded(const Ray & ray) const {
    vec3 p = ray.o - center;
    float B = dot(p, ray.d);
    float C = length2(p) - pow2(radius);
    float detSq = B * B - C;
    if ( detSq >= 0.0f ) {
        float det = std::sqrt(detSq);
        float t = - B - det;

        if ( t > ray.farT || t  < ray.nearT ) {
            t+=2 * det;
            if(t > ray.farT || t  < ray.nearT)
                return false;
        }
        return true;
    }
    return false;
}

vec3 Sphere::operator ()(Float u, Float v) const {

}

vec3 Sphere::normal(const vec3 & pos) const {
    return normalize(pos - center);
}

void Sphere::transform(const mat4 & T) {
    center = transformPoint(T, vec3(0));
    vec3  scale = extractScale(T) * vec4(vec3(1.0),0);
    radius = max(scale) * 1;
    computeBoundingBox();
    computeArea();

}

void Sphere::computeArea( ) {
    area = 4 * Constant::PI * radius * radius;
    inv_area = 1 / area;
}

void Sphere::computeBoundingBox( ) {
    BB_ = Bounds3(center - radius, center + radius);
    //todo
}

Intersection Sphere::sample(const vec2 & u, Float * pdf) const {
    Intersection it;
    it.p = center + radius * Warp::squareToUniformSphere(u);
    it.Ng=normal(it.p);
    * pdf = inv_area;
    return it;
}

Float Sphere::directPdf(const Intersection & pShape, vec3 ref) const {
    if( distance2(ref, center) <= radius * radius){
        return Primitive::directPdf(pShape,ref);
    }
    Float sinThetaMax2 = radius * radius / distance2(ref, center);
    Float cosThetaMax = std::sqrt(std::max((Float)0, 1 - sinThetaMax2));
    return 1/(2 * Constant::PI *(1-cosThetaMax));
}

Intersection Sphere::sample(const Intersection & ref, const vec2 & u, Float * pdf) const {
    //All points are visible
    Intersection it;
    if( distance2(ref.p, center) <= radius * radius){
        *pdf = 0;
        return it;
    }
    vec3 d(center-ref.p);
    Float dc = length(d);
    Float invDc = 1/dc;
    vec3 w = d * invDc;
    vec3 wx,wy;
    coordinateSystem(w,wx,wy);

    Float sinThetaMax2 = radius * radius / distance2(ref.p, center);
    Float invSinThetaMax = 1/ sqrt(sinThetaMax2);
    Float cosThetaMax = std::sqrt(std::max((Float)0, 1 - sinThetaMax2));

    Float cosTheta = (cosThetaMax-1) * u[0] +1;
    Float sinTheta2 = std::max(1-cosTheta * cosTheta,0.f);
    Float sinTheta = sqrt(sinTheta2);

    if (sinThetaMax2 < 0.00068523f /* sin^2(1.5 deg) */) {
        /* Fall back to a Taylor series expansion for small angles, where
           the standard approach suffers from severe cancellation errors */
        sinTheta2 = sinThetaMax2 * u[0];
        cosTheta = std::sqrt(1 - sinTheta2);
    }
    Float cosAlpha = sinTheta2 * invSinThetaMax +
                     cosTheta * std::sqrt(std::max((Float)0.f, 1.f - sinTheta2 * invSinThetaMax * invSinThetaMax));
    Float sinAlpha = std::sqrt(std::max((Float)0.f, 1.f - cosAlpha*cosAlpha));
    Float phi = u[1] * 2 * Constant::PI;

    vec3 nWorld = -(sinAlpha * cos(phi) * wx + sinAlpha * sin(phi) * wy +cosAlpha * w );
    vec3 pWorld = center + radius * nWorld;
    if(pWorld == ref.p){

    }

    it.p = pWorld;
    it.Ng = nWorld;
    *pdf = 1/(2 * Constant::PI *(1-cosThetaMax));
    return it;
}


