#pragma  once

#include "Integrator.hpp"
class PathIntegrator :   Integrator{
public:
    vec3 integrate(const Ray &ray, const Scene &scene, Sampler &sampler) const override;

    PathIntegrator(nlohmann::json j);
};