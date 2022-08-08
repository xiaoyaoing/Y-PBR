//2022/7/13
#pragma  once

#include "glm/glm.hpp"
#include <glm/gtx/transform.hpp>
#define  _USE_DOUBLE false
#ifdef  _USE_DOUBLE
typedef  float Float;
typedef  glm::vec2  vec2;
typedef  glm::vec3  vec3;
typedef  glm::vec4  vec4;
typedef  glm::mat4  mat4;
typedef  glm::mat3  mat3;
#else
typedef  double Float;
typedef  glm::dvec2  vec2;
typedef  glm::dvec3  vec3;
typedef  glm::dvec4  vec4;
typedef  glm::dmat3  mat3;
typedef  glm::dmat4  mat4;
#endif


typedef std::uint8_t  uint8;
typedef std::uint16_t uint16;
typedef std::uint32_t uint32;
typedef std::uint64_t uint64;

typedef std::int8_t  int8;
typedef std::int16_t int16;
typedef std::int32_t int32;
typedef std::int64_t int64;



namespace Constant{
    inline constexpr Float PI = 3.14159265358979323846;
    inline constexpr Float INV_PI = 0.31830988618379067154;
    inline constexpr Float HALF_PI = 1.57079632679489661923;
    inline constexpr Float TWO_PI = 6.283185307179586476925;
    inline constexpr Float EPSILON = 1e-9;
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









