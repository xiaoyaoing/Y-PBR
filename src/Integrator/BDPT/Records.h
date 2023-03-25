/**
 * @file 
 * @author JunPing Yuan
 * @brief 
 * @version 0.1
 * @date 2023/3/22
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma  once

#include "Common/math.hpp"

struct LightRecord{
    Float pdf;
    vec3 normal;
};

struct CameraRecord{
    vec2 pos;
};

struct SurfaceRecord{

};

struct MediumRecord{

};