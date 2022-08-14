#pragma  once

#include "Integrator.hpp"
#include "../Sampler/LightDistrib.hpp"

class PathIntegrator :   Integrator{
public:
    vec3 integrate(const Ray &ray, const Scene &scene, Sampler &sampler) const override;

    PathIntegrator(nlohmann::json j);

    virtual void  Preprocess(const Scene &scene, Sampler & sampler ) override;

private:

    std::unique_ptr<LightDistribution> lightDistribution;
    const std::string lightSampleStrategy;

};