//2022/7/13
#pragma  once

#include "Image.hpp"
#include "Ray/Ray.hpp"
#include "Common/Json.hpp"
#include "Common/Transform.hpp"

#include "Common/Json.hpp"
#include "SampleRecords/PositionAndDirectionSample.h"

class Medium;

struct PositionSample{
    vec3 p;
    Float pdf;
    Spectrum weight;
    vec3 normal;
};

struct DirSample{
    vec3 dir;
    Float pdf;
    Spectrum weight;
};

struct CameraSampleResult{

};


class Camera {
    vec3 _pos,_lookAt,_up;
    Float _fovDeg;
    Float _fovRad;
    Float _planeDist;  //distance to plane
    Float _invPlaneArea;
    Float _ratio;

    ivec2 _res;
    vec2  _pixelSize;
    //相机空间内光栅的面积
    Float A;

    mat4 _cameraToWorld;
    mat4 _toLocal;
    mat4 _rasterToCamera;
    mat4 _cameraToRaster;
public:
    Camera(const Json & c);

    void preCompute();

    PositionAndDirectionSample sampleRay(ivec2 point,vec2 posSample,vec2 dirSample) const;
    Ray sampleRay(size_t x, size_t y, /*size_t width, size_t height, Ray &ray, */vec2 sample) const;
    bool samplePosition(ivec2 point,vec2 sample,PositionSample & pSample) const;
    bool sampleDirection(ivec2 point,vec2 sample,PositionSample & pSample) const ;
    vec2 inverse(const vec3 & pos) const ;
    void drawLine(const vec3 & begin,const vec3 & end,const Spectrum & color);
    void drawLine(int x0,int x1,int y0,int y1, const Spectrum & color);
    //给定相交，采样光线，返回这条光线在胶片上的权重,光线在光栅上的位置，以及pdf
    Spectrum evalDirection(vec3 p, ivec2 * pRaster, Float * pdf) const;
    PositionAndDirectionSample sampleLi(vec3 p, ivec2 *pRaster, vec2 sample) const;
    //给定光线 返回权重 如果和光栅相交就返回光栅位置
    Spectrum rayWeight(const Ray & ray,ivec2 * pRaster) const ;
    std::shared_ptr<Medium> _medium = nullptr;
    std::unique_ptr<Image> image;

};



