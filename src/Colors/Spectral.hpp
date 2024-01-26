#pragma once

#include "Common/math.hpp"

namespace Spectral {

    static const int   CIE_samples = 471;
    static const float CIE_Min     = 360.0f;
    static const float CIE_Max     = 830.0f;

    extern const float CIE_X_entries[];
    extern const float CIE_Y_entries[];
    extern const float CIE_Z_entries[];

    void spectralXyzWeights(int samples, float lambdas[], vec3 weights[]);

    inline vec3 xyzToRgb(vec3 xyz) {
        return vec3(
            3.240479f * xyz.x + -1.537150f * xyz.y + -0.498535f * xyz.z,
            -0.969256f * xyz.x + 1.875991f * xyz.y + 0.041556f * xyz.z,
            0.055648f * xyz.x + -0.204043f * xyz.y + 1.057311f * xyz.z);
    }

    static inline vec3 wavelengthToXyz(float lambda) {
        float x    = CIE_samples * (lambda - CIE_Min) / (CIE_Max - CIE_Min);
        int   i    = clamp(int(x), 0, CIE_samples - 2);
        float u    = x - i;
        vec3  xyz0 = vec3(CIE_X_entries[i + 0], CIE_Y_entries[i + 0], CIE_Z_entries[i + 0]);
        vec3  xyz1 = vec3(CIE_X_entries[i + 1], CIE_Y_entries[i + 1], CIE_Z_entries[i + 1]);

        return xyz0 * (1.0f - u) + xyz1 * u;
    }

    static inline vec3 wavelengthToRgb(float lambda) {
        return xyzToRgb(wavelengthToXyz(lambda));
    }

}