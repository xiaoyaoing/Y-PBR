/**
 * @file 
 * @author JunPing Yuan
 * @brief 
 * @version 0.1
 * @date 2022/10/26
 *
 * @copyright Copyright (c) 2022
 *
 */
#ifndef Y_PBR_TRANSFORM_HPP
#define Y_PBR_TRANSFORM_HPP

#include "Common/math.hpp"

//use mat4 as transform

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

inline mat4 extractRotation(const mat4 & mat) {
    vec3 r = normalize(Right(mat));
    vec3 up = normalize(Up(mat));
    vec3 fwd = normalize(Forward(mat));

    return {r[0], up[0], fwd[0], 0,
            r[1], up[1], fwd[1], 0,
            r[2], up[2], fwd[2], 0,
            0, 0, 0, 1};
}

inline vec3 extractScaleVec(const mat4 & a)
{
    return vec3(length(Right(a)), length(Up(a)), length(Forward(a)));
}



////some errors in glm * ,so use this
static inline vec3 mult(const mat4 & a, const vec4 point) {
    return vec3(
            a[0][0] * point.x + a[0][1] * point.y + a[0][2] * point.z + a[0][3] * point.w,
            a[1][0] * point.x + a[1][1] * point.y + a[1][2] * point.z + a[1][3] * point.w,
            a[2][0] * point.x + a[2][1] * point.y + a[2][2] * point.z + a[2][3] * point.w
    );
}

inline vec3 transformPoint(const mat4 & a, const vec3 point) {
    return  mult(a,vec4(point,1));
}

inline vec3 transformVector(const mat4 & a, const vec3 point) {
    return mult(a,vec4(point,0));
}



inline mat4 mult(const mat4 & a, const mat4 & b) {
    mat4 result;
    for ( int i = 0 ; i < 4 ; i ++ )
        for ( int t = 0 ; t < 4 ; t ++ )
            result[i][t] =
                    a[i][0] * b[0][t] +
                    a[i][1] * b[1][t] +
                    a[i][2] * b[2][t] +
                    a[i][3] * b[3][t];

    return result;
}

inline mat4  getTranslateMatrix(const vec3 & position){
    mat4 m = mat4(1);
    m[0][3]=position.x;
    m[1][3]=position.y;
    m[2][3]=position.z;
    return m;
}

inline mat4 getScaleMatrix(const vec3 & scale){
    mat4 m = mat4(1);
    m[0][0]=scale.x;
    m[1][1]=scale.y;
    m[2][2]=scale.z;
    return m;
}

inline mat4 getTransFormMatrix(const vec3 & position, const vec3 & scale, const vec3 & rotation) {
    mat4 rotationMatrix = rotate(rotation.z, glm::dvec3(0.0, 0.0, 1.0)) *
                          rotate(rotation.y, glm::dvec3(0.0, 1.0, 0.0)) *
                          rotate(rotation.x, glm::dvec3(1.0, 0.0, 0.0));
    auto translateMatrix = getTranslateMatrix(position);
    auto scalem = getScaleMatrix(scale);
    auto  res = translateMatrix * scalem;
    return    getScaleMatrix(scale)  * rotationMatrix * translateMatrix;
}

inline mat4 getTransformNormalMat(const mat4 & matrix) {
    vec3 scaleVector(length2(Right(matrix)), length2(Up(matrix)), length2(Forward(matrix)));
    return mult(glm::scale(1.0f / scaleVector), matrix);
}

inline  mat4  getIndentifyTransform(){
    return mat4(1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                0,0,0,1);
}

namespace glm {
    void from_json(const Json & j, mat4 & transform);
}

//struct Transform {
//
//    Transform(const vec3 & position, const vec3 & scale, const vec3 & rotation);
//
//    Transform( ) : position(), scale(), rotation(), matrix(1) {
//    }
//
//
//    vec3 transformNormal(const vec3 & normal) const;
//
//    vec3 operator *(const vec3 point) const {
//        const mat4 & a = matrix;
//        return vec3(
//                a[0][0] * point.x + a[0][1] * point.y + a[0][2] * point.z + a[0][3],
//                a[1][0] * point.x + a[1][1] * point.y + a[1][2] * point.z + a[1][3],
//                a[2][0] * point.x + a[2][1] * point.y + a[2][2] * point.z + a[2][3]
//        );
//    }
//
//    vec3 transformVector(const vec3 & vec) const {
//        return mult(matrix, vec4(vec, 0));
//    }
//
//    mat4 TransformNormalMat( ) const {
//        vec3 scaleVector(length2(Right(matrix)), length2(Up(matrix)), length2(Forward(matrix)));
//        return mult(glm::scale(1.0f / scaleVector), matrix);
//    }
//
//    mat4 matrix, rotation_matrix;
//    const vec3 position, scale, rotation;
//    bool negative_determinant;
//};


#endif //Y_PBR_TRANSFORM_HPP
