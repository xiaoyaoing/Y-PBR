//2022/7/13
#pragma  once
#include "../Common/math.hpp"
#include "Image.hpp"
#include "../Ray/Ray.hpp"
#include "../Integrator/Integrator.hpp"
#include <nlohmann/json.hpp>

class Camera {
    vec3 eye;
    vec3 forward, left, up;

    Float sensor_width;
    Float focal_length;


    void samplePixel(size_t x,size_t y);

    void renderImage(Integrator * integrator,Sampler * sampler) const ;

public:
    Camera(const nlohmann::json & c);


    void sampleRay(size_t x, size_t y, size_t width, size_t height, Ray &ray, vec2 sample) const;

    void lookAt(const vec3 &p);
};



