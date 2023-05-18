#pragma  once
#include "PathVertex.h"
#include "LightPath.h"
#include "CameraPath.h"
#include "Camera/ImagePramid.h"
class BdptTracer{
  //  std::unique_ptr<PathVertex[]> lightPath;
 //   std::unique_ptr<PathVertex[]> cameraPath;
public:
    BdptTracer(const Scene & scene,const Distribution1D * distrib, const Camera * camera,ImagePramId * imagePramId = nullptr,int maxBounces=8);
    Spectrum traceSample(ivec2 pixel,Sampler & sampler);
private:
    const Scene & scene;
    const Camera * camera;
    std::unique_ptr<LightPath> cameraPath;
    std::unique_ptr<LightPath> lightPath;
    const Distribution1D * lightDistrib;
    Image * image;
    ImagePramId * imagePramid;
    int maxBounces;
};