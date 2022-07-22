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






