#pragma once
//2022/7/13
#include "../Common/math.hpp"

class Ray {

public:
    Ray();

    Ray(const Ray& rhs) = default;

    Ray(vec3 start, vec3 direction, Float nt = Constant::EPSILON, Float ft = std::numeric_limits<Float>::max());

    vec3 operator()(Float t) const {
        return o + t * d;
    }

    vec3  o;
    vec3  d;
    Float farT;
    Float nearT;
};