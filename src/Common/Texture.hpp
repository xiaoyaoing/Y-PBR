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
    virtual T eval(const Intersection * si = nullptr) const = 0;
    virtual T eval(vec2 uv) const = 0;
    virtual void makeSamplable(TextureMapJacobian jacobian) {};
    virtual vec2 sample(TextureMapJacobian jacobian, const vec2 &uv,Float * pdf) const = 0;
    virtual Float pdf(TextureMapJacobian jacobian,const vec2 & uv) const = 0;
    virtual T average() {return T();}
    virtual void setScale(Float scale){}
};

