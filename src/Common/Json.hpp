#pragma  once
#include "math.hpp"
#include "nlohmann/json.hpp"

typedef  nlohmann::json Json ;



namespace  glm{
    void from_json(const Json &j, vec3 &v);
    void from_json(const Json &j, ivec2 &v);
}


template <class T>
inline T getOptional(const Json &j, std::string field, T default_value)
{
    T ret = default_value;
    if (j.find(field) != j.end())
    {
        ret = j.at(field).get<T>();
    }
    return ret;
}

template <class T>
inline  bool containsAndGet(const Json &j, std::string field, T & value)
{
    if (j.find(field) != j.end())
    {
        value = j.at(field).get<T>();
        return true;
    }
    return false;
}

inline  bool contains(const Json &j, std::string field)
{
    return (j.find(field) != j.end());
}











