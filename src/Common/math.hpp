//2022/7/13
#pragma  once


#include "glm/glm.hpp"
#include <glm/gtx/transform.hpp>
#include "sstream"
#include "string"

#ifdef  _USE_DOUBLE

typedef  double Float;
typedef  glm::dvec2  vec2;
typedef  glm::dvec3  vec3;
typedef  glm::dvec4  vec4;
typedef  glm::dmat3  mat3;
typedef  glm::dmat4  mat4;
#else
typedef float Float;
typedef glm::vec2 vec2;
typedef glm::vec3 vec3;
typedef glm::vec4 vec4;
typedef glm::mat4 mat4;
typedef glm::mat3 mat3;
#endif

typedef glm::ivec2 ivec2;
typedef glm::ivec3 ivec3;

typedef std::uint8_t uint8;
typedef std::uint16_t uint16;
typedef std::uint32_t uint32;
typedef std::uint64_t uint64;

typedef std::int8_t int8;
typedef std::int16_t int16;
typedef std::int32_t int32;
typedef std::int64_t int64;

#define _POS_INFINY std::numeric_limits<Float>::infinity()
#define _NEG_INFINY -std::numeric_limits<Float>::infinity()
#define _NOT_IMPLEMENT_ERROR throw("This not implemented yet!");
#define _ERROR(message) throw(message);
namespace Constant {
    inline constexpr Float PI = 3.14159265358979323846;
    inline constexpr Float PiOver2 = 1.57079632679489661923;
    inline constexpr Float PiOver4 = 0.78539816339744830961;
    inline constexpr Float INV_PI = 0.31830988618379067154;
    inline constexpr Float INV_TWO_PI = INV_PI / 2;
    inline constexpr Float HALF_PI = 1.57079632679489661923;
    inline constexpr Float TWO_PI = 6.283185307179586476925;
    inline constexpr Float EPSILON = 5e-4f;
}

class   BitManip {
public :
    static inline uint32 floatBitsToUint(float f)
    {
        union {
            float f;
            uint32 i;
        } unionHack;
        unionHack.f = f;
        return unionHack.i;
    }
};

namespace std {

    template < glm::length_t Size >
    class hash < glm::vec < Size, Float, glm::defaultp>> {
    public:
        std::size_t operator ()(const glm::vec < Size, Float, glm::defaultp > & v) const {
            // See http://www.boost.org/doc/libs/1_33_1/doc/html/hash_combine.html
            uint32 result = 0;
            for ( unsigned i = 0 ; i < Size ; ++ i )
                result ^= BitManip::floatBitsToUint(v[i]) + 0x9E3779B9 + ( result << 6 ) + ( result >> 2 );
            return result;
        }
    };

}



template < glm::length_t size >
inline Float dot(glm::vec < size, Float, glm::defaultp >  a,glm::vec < size, Float, glm::defaultp > b){
    return glm::dot(a,b);
}

template < glm::length_t size >
inline  Float absDot(glm::vec < size, Float, glm::defaultp >  a,glm::vec < size, Float, glm::defaultp > b){
    return abs(dot(a,b));
}

inline Float sqr(Float a) {
    return a * a;
}

template < glm::length_t size >
inline Float length2(glm::vec < size, Float, glm::defaultp >  a) {
    return dot(a, a);
}


template < glm::length_t size >
inline Float length(glm::vec < size, Float, glm::defaultp >  a) {
    return glm::length(a);
}

template < glm::length_t size >
inline Float distance2(glm::vec < size, Float, glm::defaultp >  a,glm::vec < size, Float, glm::defaultp > b){
    return length2(a-b);
}

template < glm::length_t size >
inline Float distance(glm::vec < size, Float, glm::defaultp >  a,glm::vec < size, Float, glm::defaultp > b){
    return length(a-b);
}

inline vec3 normalize(vec3 a) {
    return glm::normalize(a);
}

inline auto radians(vec3 a) {
    return glm::radians(a);
}

inline vec3 faceForward(const vec3 & a,const vec3 & b){
    return dot(a,b)>0?a:-a;
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
   // return std::string();
    std::string s;
    std::ostringstream buffer;
    for(int i=0;i<4;i++)
        for(int j=0;j<4;j++)
            buffer<<mat[i][j]<<", ";
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

template <typename Predicate>
int FindInterval(int size, const Predicate &pred) {
    int first = 0, len = size;
    while (len > 0) {
        int half = len >> 1, middle = first + half;
        // Bisect range based on value of _pred_ at _middle_
        if (pred(middle)) {
            first = middle + 1;
            len -= half + 1;
        } else
            len = half;
    }
    return clamp(first - 1, 0, size - 2);
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
Float maxElement(const glm::vec < size, Float, glm::defaultp > & v1) {
    Float res = -1e5;
    for ( uint32 i = 0 ; i < size ; i ++ ) {
        res = std::max(v1[i], res);
    }
    return res;
}

template < glm::length_t size >
Float minElement(const glm::vec < size, Float, glm::defaultp > & v1) {
    Float res = 1e5;
    for ( uint32 i = 0 ; i < size ; i ++ ) {
        res = min(v1[i], res);
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

template < class T >
T  lerp(const T & x00, const T & x01, const T & x10, const T & x11, Float u, Float v) {

    return (x00*(1.0f - u) + x01*u)*(1.0f - v) +
           (x10*(1.0f - u) + x11*u)*v;
}

template < class T >
T  lerp(const T & x,const T & y, Float u) {

    return x * (1-u) + y * u;
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

template <class T>
T interpolate3(const T & a,const T & b,const T & c,const vec2 & uv){
    return uv.x * a + uv.y * b + (1-uv.x-uv.y)*c;
}

namespace Angle {
    static Float radToDeg(Float a) {
        return a * ( 180.0 / Constant::PI );
    }

    static Float degToRad(Float a) {
        return a * ( Constant::PI / 180.0 );
    }
}



inline Float PowerHeuristic(Float a,Float b){
    return a*a/(a*a + b*b);
}











