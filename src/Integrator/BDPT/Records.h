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
#include "Camera/Camera.hpp"
#include "Ray/Intersection.hpp"
#include "Lights/Light.hpp"
#include "SampleRecords/PositionAndDirectionSample.h"

struct LightRecord{
    PositionAndDirectionSample sample;
    float lightPdf;
};

struct CameraRecord{
    ivec2 pixel;
    PositionAndDirectionSample sample;
};

struct SurfaceRecord{
    SurfaceEvent event;
};

struct MediumRecord{

};