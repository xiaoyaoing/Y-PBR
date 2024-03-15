#include "../Common/math.hpp"
#pragma once
using Spectrum = vec3;

inline bool isBlack(const Spectrum& color) {
    static const Float eps = 1e-5f;
    return color.x <= eps && color.y <= eps && color.z <= eps;
}

inline Float luminace(const Spectrum& color) {
    const Float YWeight[3] = {0.212671f, 0.715160f, 0.072169f};
    return YWeight[0] * color[0] + YWeight[1] * color[1] + YWeight[2] * color[2];
}

bool hasNan(const Spectrum& color);

inline bool hasNeg(const Spectrum& color) {
    return color[0] < 0 || color[1] < 0 || color[2] < 0;
}

//void from_json(const Json &j, Spectrum & spectrum);