#pragma once

#include "Integrator/PathIntegrator.hpp"
#include "work-queue.hpp"

struct Bucket
    {
        Bucket() : min(0), max(0) { }
        Bucket(const ivec2& min, const ivec2& max) : min(min), max(max) { }

        ivec2 min;
        ivec2 max;
    };


class Render {
public:
    Render(const Json & json);
    void Go();
private:
    std::shared_ptr<Camera> camera ;
    std::shared_ptr<Scene>  scene;
    std::unique_ptr<Integrator> integrator;
    std::shared_ptr<Sampler> sampler;
};