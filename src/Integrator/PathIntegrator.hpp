#pragma  once

#include "../Sampler/LightDistrib.hpp"
#include "AbstractPathTracer.hpp"
struct PathTraceSettings {
    bool enableTwoSide;
};

class PathIntegrator :   AbstractPathTracer{
public:
    vec3 integrate(const Ray &ray, const Scene &scene, Sampler &sampler) const override;

    PathIntegrator(nlohmann::json j);

    virtual void  Preprocess(const Scene &scene, Sampler & sampler ) override;

private:

    std::unique_ptr<LightDistribution> lightDistribution;
    const std::string lightSampleStrategy;

};