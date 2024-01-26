#pragma once

#include "scene.hpp"
#include "Ray/Ray.hpp"
#include "Camera/Camera.hpp"
#include "Sampler/Distrib.hpp"
#include "Sampler/Sampler.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"

class Medium;

//class Sampler;

//class Image;

class Integrator {
public:
    Integrator(const Json& json) {
        minBounces = getOptional(json, "min_bounces", 0);
        maxBounces = getOptional(json, "max_bounces", 8);
    }
    Integrator() : minBounces(0), maxBounces(8) {
    }
    virtual void render(const Scene& scene) = 0;

    virtual void process(const Scene& scene, Sampler& sampler) = 0;

protected:
    int minBounces, maxBounces;
};

//Integrator Based Sampler
class SamplerIntegrator : public Integrator {

public:
    SamplerIntegrator(std::shared_ptr<Camera> camera, std::shared_ptr<Sampler> sampler, const Json& json) : _camera(camera), _sampler(sampler), Integrator(json) {}
    SamplerIntegrator(std::shared_ptr<Camera> camera, std::shared_ptr<Sampler> sampler) : _camera(camera), _sampler(sampler) {}

    void         render(const Scene& scene) override;
    virtual vec3 integrate(const Ray& ray, const Scene& scene, Sampler& sampler) const = 0;
    virtual void renderPixel(int x, int y) const;

public:
    std::shared_ptr<Camera>  _camera;
    std::shared_ptr<Sampler> _sampler;
    int                      tileISze = 16;
};