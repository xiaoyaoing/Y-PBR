#include "Integrator.hpp"
#include "PathIntegrator.hpp"

class VolPathIntegrator : public PathIntegrator {
public:
    VolPathIntegrator(std::shared_ptr<Camera> camera, std::shared_ptr<Sampler> sampler, const Json& json) : PathIntegrator(camera, sampler, json) {}

    vec3 integrate(const Ray& ray, const Scene& scene, Sampler& sampler) const override;
};