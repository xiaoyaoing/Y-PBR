#pragma once
#include "Common/math.hpp"
#include "Ray/Intersection.hpp"

class SurfaceIntersection;

enum TextureMapJacobian {
    MAP_UNIFORM,
    MAP_SPHERICAL,
    MAP_JACOBIAN_COUNT,
};

template<class T>
class Texture{
public:
    virtual T Evaluate(const Intersection * si = nullptr) const = 0;
    virtual T Evaluate(const vec2 & uv) const = 0;
    virtual void makeSamplable(TextureMapJacobian jacobian) {};
    virtual vec2 sample(TextureMapJacobian jacobian, const vec2 &uv) const = 0;
    virtual Float pdf(TextureMapJacobian jacobian,const vec2 & uv) const = 0;
};

