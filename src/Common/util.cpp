#include "util.hpp"

void glm::from_json(const nlohmann::json &j, vec3 & v)
{
    if (j.type() == nlohmann::json::value_t::array)
        for (int i = 0; i < 3; i++) j.at(i).get_to(v[i]);
    else
        for (int i = 0; i < 3; i++) j.get_to(v[i]);
}

Transform::Transform(const vec3 &position, const vec3 &scale, const vec3 &rotation)
    :position(position),scale(scale),rotation(rotation)
{
    rotation_matrix = rotate(rotation.z, glm::dvec3(0.0, 0.0, 1.0)) *
                      rotate(rotation.y, glm::dvec3(0.0, 1.0, 0.0)) *
                      rotate(rotation.x, glm::dvec3(1.0, 0.0, 0.0));

    matrix = glm::translate(mat4(1.0), position) *
             rotation_matrix *
             glm::scale(mat4(1.0), scale);

}

vec3 Transform::transformNormal(const vec3 &normal) const {

}
