#include "Json.hpp"
#include "Common/Transform.hpp"

void glm::from_json(const Json& j, vec3& v) {
    if (j.type() == Json::value_t::array)
        for (int i = 0; i < 3; i++) j.at(i).get_to(v[i]);
    else
        for (int i = 0; i < 3; i++) j.get_to(v[i]);
}

void glm::from_json(const Json& j, ivec2& v) {
    if (j.type() == Json::value_t::array)
        for (int i = 0; i < 2; i++) j.at(i).get_to(v[i]);
    else
        for (int i = 0; i < 2; i++) j.get_to(v[i]);
}