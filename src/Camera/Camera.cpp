//2022/7/13

#include "Camera.hpp"
#include "../Common/util.hpp"
Camera::Camera(const nlohmann::json & c) {
    eye = c.at("eye");
//    from_json(c.at("eye"),eye);
    focal_length = c.at("focal_length").get<Float>() / Float(1000);
    sensor_width = c.at("sensor_width").get<Float>() / Float(1000);
    sample_count = c.at("sample_count").get<uint32>();

    vec3 look_at=c.at("look_at");
    lookAt(look_at);
}

void Camera::sampleRay(size_t x, size_t y ,size_t width,size_t height, Ray & ray, vec2 sample) const {

    Float pixel_size = sensor_width / width;

    vec2 half_dim = vec2(width, height) * Float(0.5);

    vec2 px(x + sample[0], y + sample[1]);
    vec2 local = pixel_size *  (half_dim - px);
    vec3 direction = glm::normalize(forward * focal_length + left * local.x + up * local.y);
    // Pinhole camera ray
    ray= Ray(eye, direction);
}



void Camera::samplePixel(size_t x, size_t y) {

}

void Camera::renderImage(Integrator *integrator, Sampler *sampler) const {

}

void Camera::lookAt(const vec3& p)
{
    forward = normalize(p - eye);
    left = cross({0.0, 1.0, 0.0}, forward);
    left = length(left) < Constant::EPSILON ? vec3(-1.0, 0.0, 0.0) : normalize(left);
    up = normalize(glm::cross(forward, left));
}








