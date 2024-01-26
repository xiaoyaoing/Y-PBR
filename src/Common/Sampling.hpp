#pragma once

namespace Warper {

    vec3 UniformSampleSphere(const Point2f& u) {
        Float z   = 1 - 2 * u[0];
        Float r   = std::sqrt(std::max((Float)0, (Float)1 - z * z));
        Float phi = 2 * Pi * u[1];
        return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
    }
}