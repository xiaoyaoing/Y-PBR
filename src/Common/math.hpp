//2022/7/13
#pragma  once

#include "glm/glm.hpp"
#include <glm/gtx/transform.hpp>
#define  _USE_DOUBLE false
#ifdef  _USE_DOUBLE
typedef  float Float;
typedef  glm::vec2  vec2;
typedef  glm::vec3  vec3;
typedef  glm::mat4  mat4;
typedef  glm::mat3  mat3;
#else
typedef  double Float;
typedef  glm::dvec2  vec2;
typedef  glm::dvec3  vec3;
typedef  glm::dmat3  mat3;
typedef  glm::dmat4  mat4;
#endif

namespace Constant{
    inline constexpr double PI = 3.14159265358979323846;
    inline constexpr double INV_PI = 0.31830988618379067154;
    inline constexpr double HALF_PI = 1.57079632679489661923;
    inline constexpr double TWO_PI = 6.283185307179586476925;
    inline constexpr double EPSILON = 1e-9;
}

inline  Float dot(vec3 a,vec3 b){
    return glm::dot(a,b);
}

inline  Float pow2(Float a){
    return a*a;
}

inline  Float length2(vec3 a){
    return dot(a,a);
}

inline  vec3 normalize(vec3 a){
    return glm::normalize(a);
}

inline  auto radians(vec3 a){
    return glm::radians(a);
}

inline mat4 rotate(Float radians,vec3 axis){
    return glm::rotate(radians,axis);
}


inline bool solveQuadratic(double a, double b, double c, double& t_min, double& t_max)
{
    if (a != 0.0)
    {
        double d = b * b - 4.0 * a * c;

        if (d < 0.0) return false;

        double t = -0.5 * (b + (b < 0.0 ? -std::sqrt(d) : std::sqrt(d)));

        t_min = t / a;
        t_max = c / t;

        if (t_min > t_max) std::swap(t_min, t_max);

        return true;
    }
    if (b != 0.0)
    {
        t_min = t_max = -c / b;
        return true;
    }
    return false;
}




