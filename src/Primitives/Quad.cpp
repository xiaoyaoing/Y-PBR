//2022/7/17

#include "Quad.hpp"
#include "Bsdfs/Reflection.hpp"
#include "iostream"
#include "Common/Transform.hpp"


void Quad::preCompute( ) {
    computeArea();
    computeBoundingBox();
    _invSq = vec2(1 / length2(_edge0), 1 / length2(_edge1));
    _normal = normalize(cross(_edge1, _edge0));
}


std::optional < Intersection > Quad::intersect(Ray & ray) const {
    Intersection its;
    vec3 n = normal(_base);
    Float dotW = dot(ray.d, n);
    if ( std::abs(dotW) < Constant::EPSILON ) {
        return std::nullopt;
    }
    Float t = dot(n, ( _base - ray.o )) / dotW;

    if ( t < ray.nearT || t > ray.farT )
        return std::nullopt;

    vec3 q = ray(t);
    vec3 v = q - _base;
    Float l0 = dot(v, _edge0) * _invSq.x;
    Float l1 = dot(v, _edge1) * _invSq.y;


    if ( l0 < 0.0f || l0 > 1.0f || l1 < 0.0f || l1 > 1.0f )
        return std::nullopt;

    ray.farT = t;
    its.uv = vec2(l0, l1);
    its.p = q;
   its.Ng = its.Ns =  _normal;
    its.primitive = this;
    its.bsdf = bsdf.get();
    return {its};
}

bool Quad::occluded(const Ray & ray) const {
    vec3 n = normal(_base);
    Float dotW = dot(ray.d, n);
    if ( std::abs(dotW) < Constant::EPSILON ) {
        return false;
    }
    Float t = dot(n, ( _base - ray.o )) / dotW;
    if ( t < ray.nearT || t > ray.farT )
        return false;

    vec3 q = ray(t);
    vec3 v = q - _base;
    Float l0 = dot(v, _edge0) * _invSq.x;
    Float l1 = dot(v, _edge1) * _invSq.y;


    if ( l0 < 0.0f || l0 > 1.0f || l1 < 0.0f || l1 > 1.0f )
        return false;
    return true;
}


Intersection Quad::sample(const vec2 & u, Float * pdf) const {
    Intersection its;
    its.p = _base + _edge0 * u[0] + _edge1 * u[1];
    its.Ng = normal(its.p);
    * pdf = inv_area;
    return its;
}


vec3 Quad::normal(const vec3 & pos) const {
    return _normal;
}

void Quad::transform(const mat4 & T) {

    _base = transformPoint(T,vec3(0.0f));
    _edge0 = transformVector(T,vec3(1.0f, 0.0f, 0.0f));
    _edge1 = transformVector(T,vec3(0.0f, 0.0f, 1.0f));

    _base -= _edge0 * 0.5f;
    _base -= _edge1 * 0.5f;

    preCompute();
}

void Quad::computeBoundingBox( ) {

    Bounds3 bounds;
    bounds = Union(bounds, _base);
    bounds = Union(bounds, _base + _edge1);
    bounds = Union(bounds, _base + _edge0);
    bounds = Union(bounds, _base + _edge0 + _edge1);
    BB_ = bounds;
}

void Quad::computeArea( ) {
    area = glm::length(cross(_edge0, _edge1));
    inv_area = 1 / area;
}

Float Quad::directPdf(const Intersection & pShape, vec3 ref) const {
    return length2(ref-pShape.p) * InvArea() / dot(-pShape.w,pShape.Ng) ;
}

Float Quad::powerToRadianceScale( ) const {
    return Constant::INV_PI * inv_area;
}


