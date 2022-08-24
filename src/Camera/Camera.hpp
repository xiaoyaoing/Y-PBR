//2022/7/13
#pragma  once
#include "../Common/math.hpp"
#include "Image.hpp"
#include "../Ray/Ray.hpp"
#include "../Integrator/Integrator.hpp"
#include <nlohmann/json.hpp>

class Camera {
    vec3 _pos,_lookAt,_up;
    Float _fovDeg;
    Float _fovRad;
    Float _planeDist;  //distance to plane
    Float _invPlaneArea;
    Float _ratio;

    ivec2 _res;
    vec2  _pixelSize;

    std::unique_ptr<Image> image;

    Transform _transform;


    void samplePixel(size_t x,size_t y);

    void renderImage(Integrator * integrator,Sampler * sampler) const ;

public:
    Camera(const nlohmann::json & c);

    void preCompute();

    void sampleRay(size_t x, size_t y, /*size_t width, size_t height,*/ Ray &ray, vec2 sample) const;

   // void lookAt(const vec3 &p);

    uint32 sample_count  ;

    void Capture();


};



