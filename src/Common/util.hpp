#pragma  once
#include "math.hpp"
#include "nlohmann/json.hpp"


namespace  glm{
    void from_json(const nlohmann::json &j, vec3 &v);
    void from_json(const nlohmann::json &j, ivec2 &v);
}



struct Transform
{
    Transform(const vec3 &position, const vec3 &scale, const vec3 &rotation);

    Transform():position(),scale(),rotation(),matrix(1) {
    }


    vec3 transformNormal(const vec3 & normal) const;

    vec3 operator * (const vec3 point) const {
        const mat4 & a =matrix;
        return vec3(
                a[0][0]*point.x + a[0][1]*point.y + a[0][2]*point.z + a[0][3],
                a[1][0]*point.x + a[1][1]*point.y + a[1][2]*point.z + a[1][3],
                a[2][0]*point.x + a[2][1]*point.y + a[2][2]*point.z + a[2][3]
        );
    }

    vec3 transformVector(const vec3 & vec) const {
        return mult(matrix,vec4(vec,0));
    }

    mat4 TransformNormalMat() const {
        vec3 scaleVector(length2(Right(matrix)), length2(Up(matrix)), length2(Forward(matrix)));
        return mult(glm::scale(1.0f/scaleVector),matrix);
    }

    mat4 matrix, rotation_matrix;
    const vec3  position, scale, rotation;
    bool negative_determinant;
};

void from_json(const nlohmann::json & j,Transform & transform);


template <class T>
inline T getOptional(const nlohmann::json &j, std::string field, T default_value)
{
    T ret = default_value;
    if (j.find(field) != j.end())
    {
        ret = j.at(field).get<T>();
    }
    return ret;
}

template <class T>
inline  bool containsAndGet(const nlohmann::json &j, std::string field, T & value)
{
    if (j.find(field) != j.end())
    {
        value = j.at(field).get<T>();
        return true;
    }
    return false;
}

inline  bool contains(const nlohmann::json &j, std::string field)
{
    return (j.find(field) != j.end());
}





/// Simple floating point clamping function


/// Simple integer clamping function
inline int clamp(int value, int min, int max) {
    if (value < min)
        return min;
    else if (value > max)
        return max;
    else return value;
}








