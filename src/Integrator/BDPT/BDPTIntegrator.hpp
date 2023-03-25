#pragma once

#include "Integrator/Integrator.hpp"
#include "PathVertex.h"
int generateLightPath(const Distribution1D * ptr, Sampler & sampler,int maxDepth,PathVertex * path);
int generateCameraPath();
int randomWalk(const Scene &scene, Ray &ray, Sampler &sampler, int maxDepth, PathVertex *pVertex);
void connectPath();
class BDPTIntegrator : SamplerIntegrator{
    void process(const Scene &scene, Sampler &sampler) override;

    vec3 integrate(const Ray &ray, const Scene &scene, Sampler &sampler) const override;

    void render(const Scene &scene) const override;

protected:
    std::unique_ptr<Distribution1D> lightDistrib;
};