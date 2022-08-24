//2022/7/13
#pragma  once


#include "glm/glm.hpp"
#include <glm/gtx/transform.hpp>
#include "sstream"
#include "string"

#define  _USE_DOUBLE false
#ifdef  _USE_DOUBLE
typedef float Float;
typedef glm::vec2 vec2;
typedef glm::vec3 vec3;
typedef glm::vec4 vec4;
typedef glm::mat4 mat4;
typedef glm::mat3 mat3;
typedef glm::ivec2 ivec2;
#else
typedef  double Float;
typedef  glm::dvec2  vec2;
typedef  glm::dvec3  vec3;
typedef  glm::dvec4  vec4;
typedef  glm::dmat3  mat3;
typedef  glm::dmat4  mat4;
#endif

typedef glm::ivec2 ivec2;

typedef std::uint8_t uint8;
typedef std::uint16_t uint16;
typedef std::uint32_t uint32;
typedef std::uint64_t uint64;

typedef std::int8_t int8;
typedef std::int16_t int16;
typedef std::int32_t int32;
typedef std::int64_t int64;


namespace Constant {
    inline constexpr Float PI = 3.14159265358979323846;
    inline constexpr Float INV_PI = 0.31830988618379067154;
    inline constexpr Float HALF_PI = 1.57079632679489661923;
    inline constexpr Float TWO_PI = 6.283185307179586476925;
    inline constexpr Float EPSILON = 1e-9;
}

inline Float dot(vec3 a, vec3 b) {
    return glm::dot(a, b);
}

inline Float pow2(Float a) {
    return a * a;
}

inline Float length2(vec3 a) {
    return dot(a, a);
}


inline vec3 normalize(vec3 a) {
    return glm::normalize(a);
}

inline auto radians(vec3 a) {
    return glm::radians(a);
}

inline vec3 intToColor(uint32_t i) {
    return vec3(( i >> 16 ) & 0xFF, ( i >> 8 ) & 0xFF, i & 0xFF) / (Float) 255.0;
}

inline mat4 rotate(Float radians, vec3 axis) {
    return glm::rotate(radians, axis);
}

inline Float clamp(Float value, Float min, Float max) {
    if ( value < min )
        return min;
    else if ( value > max )
        return max;
    else return value;
}

inline std::string toColorStr(const vec3 & v) {
    std::string s;
    std::ostringstream buffer;
    buffer << "r:" << v.x << " g: " << v.y << " b: " << v.z;
    return buffer.str();
}

inline std::string Mat4ToStr(const mat4 & mat) {
    std::string s;
    std::ostringstream buffer;
    for(int i=0;i<4;i++)
        for(int j=0;j<4;j++)
            buffer<<mat[i][j]<<" ";
    buffer<<std::endl;
    return buffer.str();
}


template < int num >
Float max(const glm::vec < num, Float, glm::defaultp > & v) {
    Float maxVal = -1e5;
    for ( uint32 i = 0 ; i < num ; i ++ ) {
       maxVal = std::max(v[i],maxVal);
    }
    return maxVal;
}


template < uint32 num >
bool AllInRange(const glm::vec < num, Float, glm::defaultp > & v,
                Float minVal, Float maxVal) {
    for ( uint32 i = 0 ; i < num ; i ++ ) {
        if ( v[i] < minVal || v[i] > maxVal ) {
            return false;
        }
    }

    return true;
}

template < glm::length_t size >
glm::vec < size, Float, glm::defaultp > pow(const glm::vec < size, Float, glm::defaultp > & v,
                                            Float e) {
    glm::vec < size, Float, glm::defaultp > res;
    for ( uint32 i = 0 ; i < size ; i ++ ) {
        res[i] = std::pow(v[i], e);
    }
    return res;
}


template < glm::length_t size >
glm::vec < size, Float, glm::defaultp > max(const glm::vec < size, Float, glm::defaultp > & v1,
                                            const glm::vec < size, Float, glm::defaultp > & v2) {
    glm::vec < size, Float, glm::defaultp > res;
    for ( uint32 i = 0 ; i < size ; i ++ ) {
        res[i] = v1[i] > v2[i] ? v1[i] : v2[i];
    }
    return res;
}

template < glm::length_t size >
glm::vec < size, Float, glm::defaultp > min(const glm::vec < size, Float, glm::defaultp > & v1,
                                            const glm::vec < size, Float, glm::defaultp > & v2) {
    glm::vec < size, Float, glm::defaultp > res;
    for ( uint32 i = 0 ; i < size ; i ++ ) {
        res[i] = v1[i] < v2[i] ? v1[i] : v2[i];
    }
    return res;
}

template < glm::length_t size >
glm::vec < size, Float, glm::defaultp > clamp(const glm::vec < size, Float, glm::defaultp > & v, Float l, Float r) {
    glm::vec < size, Float, glm::defaultp > res;
    for ( uint32 i = 0 ; i < size ; i ++ ) {
        res[i] = clamp(v[i], l, r);
    }
    return res;
}

template < glm::length_t size >
glm::vec < size, Float, glm::defaultp > clamp(const glm::vec < size, Float, glm::defaultp > & v,
                                              const glm::vec < size, Float, glm::defaultp > & l,
                                              const glm::vec < size, Float, glm::defaultp > & r
                                              ) {
    glm::vec < size, Float, glm::defaultp > res;
    for ( uint32 i = 0 ; i < size ; i ++ ) {
        res[i] = clamp(v[i], l[i], r[i]);
    }
    return res;
}

template < glm::length_t size >
glm::vec < size, Float, glm::defaultp > lerp(const glm::vec < size, Float, glm::defaultp > & v1,
                                             const glm::vec < size, Float, glm::defaultp > & v2,
                                             const glm::vec < size, Float, glm::defaultp > & p
) {
    glm::vec < size, Float, glm::defaultp > res;
    for ( uint32 i = 0 ; i < size ; i ++ ) {
        res[i] = std::lerp(v1[i], v2[i], p);
    }
    return res;
}

template < glm::length_t size >
Float compAdd(const glm::vec < size, Float, glm::defaultp > & v) {
    Float res = 0;
    for ( uint32 i = 0 ; i < size ; i ++ ) {
        res += v[i];
    }
    return res;
}

template < glm::length_t size >
bool hasNav(const glm::vec < size, Float, glm::defaultp > & v) {
    Float res;
    for ( uint32 i = 0 ; i < size ; i ++ ) {
        if ( v[i] < 0 )
            return true;
    }
    return false;
}

template < glm::length_t size >
int maxDim(const glm::vec < size, Float, glm::defaultp > & v) {
    Float maxVal = -1e5;
    int dim=-1;
    for ( uint32 i = 0 ; i < size ; i ++ ) {
        if(v[i]>maxVal){
            maxVal = v[i];
            dim=int(i);
        }
    }
    return dim;
}

namespace Angle {
    static Float radToDeg(Float a) {
        return a * ( 180.0 / Constant::PI );
    }

    static Float degToRad(Float a) {
        return a * ( Constant::PI / 180.0 );
    }
}

inline vec3 Right(const mat4 & m) {
    return vec3(m[0][0], m[1][0], m[2][0]);
}

inline vec3 Up(const mat4 & m) {
    return vec3(m[0][1], m[1][1], m[2][1]);
}

inline vec3 Forward(const mat4 & m) {
    return vec3(m[0][2], m[1][2], m[2][2]);
}

inline mat4 extractScale(const mat4 & mat) {
    return glm::scale(vec3(glm::length(Right(mat)), glm::length(Up(mat)), glm::length(Forward(mat))));
}

inline mat4 extractRotation(const mat4 & mat){
    vec3 r = normalize(Right(mat));
    vec3 up = normalize(Up(mat));
    vec3 fwd= normalize(Forward(mat));

    return {r[0],up[0],fwd[0],0,
                r[1],up[1],fwd[1],0,
                r[2],up[2],fwd[2],0,
                0,0,0,1};
}


//some errors in glm * ,so use this
inline vec3 mult(const mat4 & a ,const vec4  point ){
    return vec3(
            a[0][0]*point.x + a[0][1]*point.y + a[0][2]*point.z + a[0][3] * point.w,
            a[1][0]*point.x + a[1][1]*point.y + a[1][2]*point.z + a[1][3] * point.w,
            a[2][0]*point.x + a[2][1]*point.y + a[2][2]*point.z + a[2][3] * point.w
    );
}











