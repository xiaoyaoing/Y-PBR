#pragma once
#include "PathVertex.h"
#include "LightPath.h"
#include "Camera/ImagePrid.h"
class BdptTracer {
    //  std::unique_ptr<PathVertex[]> lightPath;
    //   std::unique_ptr<PathVertex[]> cameraPath;
public:
    BdptTracer(const Scene& scene, const Distribution1D& distrib, const std::map<const Light*, size_t>& lightIdxs, const Camera* camera, ImagePrid* imagePramId = nullptr, int maxBounces = 8);
    Spectrum traceSample(ivec2 pixel, Sampler& sampler);

private:
    const Scene&                          scene;
    const Camera*                         camera;
    std::unique_ptr<LightPath>            cameraPath;
    std::unique_ptr<LightPath>            lightPath;
    const Distribution1D&                 lightDistrib;
    const std::map<const Light*, size_t>& lightIdxs;
    Image*                                image;
    ImagePrid*                            imagePramid;
    int                                   maxBounces;
};