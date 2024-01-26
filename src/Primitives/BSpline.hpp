#pragma once
#include "Common/math.hpp"

namespace BSpline {

    // http://www.answers.com/topic/b-spline
    template<typename T>
    inline T quadratic(T p0, T p1, T p2, Float t) {
        return (0.5f * p0 - p1 + 0.5f * p2) * t * t + (p1 - p0) * t + 0.5f * (p0 + p1);
    }
    template<typename T>
    inline T quadratic(T p0, T p1, T p2, T p3, Float t) {
        return quadratic(lerp(p0, p1, t), lerp(p1, p2, t), lerp(p2, p3, t), t);
    }

    template<typename T>
    inline T quadraticDeriv(T p0, T p1, T p2, Float t) {
        return (p0 - 2.0f * p1 + p2) * t + (p1 - p0);
    }

    template<typename T>
    inline T quadraticDeriv(T p0, T p1, T p2, T p3, Float t) {
        T p00 = p1 - p0;
        T p01 = p2 - p1;
        T p02 = p3 - p2;

        const float t0 = 1.0f - t, t1 = t;
        const float B0 = -(t0 * t0);
        const float B1 = -2.0f * (t0 * t1) + (t0 * t0);
        const float B2 = +2.0f * (t0 * t1) - (t1 * t1);
        const float B3 = +(t1 * t1);
        return 3.0f * (B0 * p0 + B1 * p1 + B2 * p2 + B3 * p3);

        return quadraticDeriv(p00, p01, p02, t);
    }

    inline vec2 quadraticMinMax(Float p0, Float p1, Float p2) {
        Float xMin = (p0 + p1) * 0.5f;
        Float xMax = (p1 + p2) * 0.5f;
        if (xMin > xMax)
            std::swap(xMin, xMax);

        Float tFlat = (p0 - p1) / (p0 - 2.0f * p1 + p2);
        if (tFlat > 0.0f && tFlat < 1.0f) {
            Float xFlat = quadratic(p0, p1, p2, tFlat);
            xMin        = std::min(xMin, xFlat);
            xMax        = std::max(xMax, xFlat);
        }
        return vec2(xMin, xMax);
    }

}