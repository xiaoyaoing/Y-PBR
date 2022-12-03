//2022/7/13
#pragma  once

#include "Image.hpp"
#include "Ray/Ray.hpp"
#include "Common/Json.hpp"
#include "Common/Transform.hpp"

#include "Common/Json.hpp"

class Medium;

class Camera {
    vec3 _pos,_lookAt,_up;
    Float _fovDeg;
    Float _fovRad;
    Float _planeDist;  //distance to plane
    Float _invPlaneArea;
    Float _ratio;

    ivec2 _res;
    vec2  _pixelSize;

    mat4 _toWorld;
    mat4 _toLocal;
public:
    Camera(const Json & c);

    void preCompute();

    Ray sampleRay(size_t x, size_t y, /*size_t width, size_t height, Ray &ray, */vec2 sample) const;
    vec2 inverse(const vec3 & pos) const ;
    void drawLine(const vec3 & begin,const vec3 & end,const Spectrum & color);
    void drawLine(int x0,int x1,int y0,int y1, const Spectrum & color);

    std::shared_ptr<Medium> _medium = nullptr;
    std::unique_ptr<Image> image;
};



