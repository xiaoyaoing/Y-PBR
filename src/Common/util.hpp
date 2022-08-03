#pragma  once
#include "math.hpp"
#include "nlohmann/json.hpp"



namespace  glm{
    void from_json(const nlohmann::json &j, vec3 &v);
}


struct Transform
{
    Transform(const vec3 &position, const vec3 &scale, const vec3 &rotation);

    vec3 transformNormal(const vec3 & normal) const;

    vec3 operator()(const vec3 point) const {
        return matrix * vec4(point,1.0);
    }

    mat4 matrix, rotation_matrix;
    const vec3  position, scale, rotation;
    bool negative_determinant;
};

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

/// Simple floating point clamping function
inline float clamp(Float value, Float min, Float max) {
    if (value < min)
        return min;
    else if (value > max)
        return max;
    else return value;
}

/// Simple integer clamping function
inline int clamp(int value, int min, int max) {
    if (value < min)
        return min;
    else if (value > max)
        return max;
    else return value;
}








