#pragma  once

#include "Sampler/LightDistrib.hpp"
#include "Integrator.hpp"

struct PathTraceSettings {
    bool enableTwoSide;
};

class PathIntegrator :  public SamplerIntegrator{
public:
    PathIntegrator(std::shared_ptr<Camera> camera,std::shared_ptr<Sampler> sampler) : SamplerIntegrator(camera,sampler){}
    vec3 integrate(const Ray &ray, const Scene &scene, Sampler &sampler) const override;
//    PathIntegrator(Json j);
    virtual void  process(const Scene &scene, Sampler & sampler ) override;

protected:
    std::unique_ptr<LightDistribution> lightDistribution;
    const std::string lightSampleStrategy;
};