#include "Integrator.hpp"
class VolPathIntegrator : SamplerIntegrator{
    void process(const Scene & scene, Sampler & sampler) override;

    vec3 integrate(const Ray & ray, const Scene & scene, Sampler & sampler) const override;
protected :
    std::unique_ptr<Distribution1D> lightDistribution;
};