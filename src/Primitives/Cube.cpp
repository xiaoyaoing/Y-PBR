#include "Cube.hpp"
#include "spdlog/spdlog.h"
#include "Common/Transform.hpp"



std::optional < Intersection > Cube::intersect(Ray & ray) const {
    Intersection its;

    vec3 p = mult(_invRot ,vec4((ray.o - _pos),1));
    vec3 d = mult(_invRot , vec4(ray.d,1));

    vec3 invD = 1.0f /  d;
    vec3 relMin((-_scale - p));
    vec3 relMax(( _scale - p));

    float ttMin = ray.nearT, ttMax = ray.farT;
    for (int i = 0; i < 3; ++i) {
        if (invD[i] >= 0.0f) {
            ttMin = std::max(ttMin, relMin[i]*invD[i]);
            ttMax = std::min(ttMax, relMax[i]*invD[i]);
        } else {
            ttMax = std::min(ttMax, relMin[i]*invD[i]);
            ttMin = std::max(ttMin, relMax[i]*invD[i]);
        }
    }

    if (ttMin <= ttMax) {
        Float t=1.0 ;
        if (ttMin > ray.nearT && ttMin < ray.farT) {
            t = ttMin;
        }
        else  if(ttMax > ray.nearT && ttMax < ray.farT)
        {
            t = ttMax;
        }
        ray.farT = t;
        its.p= ray(t);
        its.Ns = its.Ng = normal(its.p);
        its.primitive = this;
        its.bsdf=bsdf.get();
        return {its};
    }
    return std::nullopt;
}

bool Cube::occluded(const Ray & ray) const {

    vec3 p = transformVector(_invRot ,ray.o - _pos);
    vec3 d = transformVector(_invRot , ray.d);

    vec3 invD = 1.0f /  d;
    vec3 relMin((-_scale - p));
    vec3 relMax(( _scale - p));

    float ttMin = ray.nearT, ttMax = ray.farT;
    for (int i = 0; i < 3; ++i) {
        if (invD[i] >= 0.0f) {
            ttMin = std::max(ttMin, relMin[i]*invD[i]);
            ttMax = std::min(ttMax, relMax[i]*invD[i]);
        } else {
            ttMax = std::min(ttMax, relMin[i]*invD[i]);
            ttMin = std::max(ttMin, relMax[i]*invD[i]);
        }
    }
    return ttMin <= ttMax;
}

vec3 Cube::normal(const vec3 & pos) const {
    vec3 p = transformPoint(_invRot ,pos-_pos);
    vec3  n(0.0f);
    int dim = maxDim((abs(p) - _scale));
    n[dim] = p[dim] < 0.0f ? -1.0f : 1.0f;
    n = mult(_rot , vec4(n,1));
    return n;
}

void Cube::computeArea( ) {
    area = 8 *(_scale.x*_scale.y+_scale.y*_scale.z+_scale.z*_scale.x);
    inv_area=1/area;
}

void Cube::computeBoundingBox( ) {
    Bounds3 box;
    for (int i = 0; i < 8; ++i) {
        vec3 v= mult( _rot ,vec4((i & 1 ? _scale.x : -_scale.x),
                            (i & 2 ? _scale.y : -_scale.y),
                            (i & 4 ? _scale.z : -_scale.z),1));
        box= Union(box,_pos + v);;
    }

    BB_ = box;
}

void Cube::transform(const mat4 & T) {
    _pos = transformPoint(T,vec3(0.0f));
    _scale = extractScale(T) * vec4(vec3(0.5f),1);
    _rot=extractRotation(T);
    _invRot = glm::transpose(_rot);
    auto s1 = Mat4ToStr(_rot);
    auto s2 = Mat4ToStr(_invRot);
    preCompute();
}

vec3 Cube::operator ()(Float u, Float v) const {
    return vec3(0);
}

Intersection Cube::sample(const vec2 & u, Float * pdf) const {
    return Intersection();
    //to do
}