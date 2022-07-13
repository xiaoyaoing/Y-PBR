//2022/7/13
#include "glm/glm.hpp"
#include <nlohmann/json.hpp>

class Camera {
    glm::dvec3 eye;
    glm::dvec3 forward, left, up;

    Camera(const nlohmann::json &j);
};



