#pragma once

#include "Integrator/Integrator.hpp"
#include "BdptTracer.h"

///Bidrectional path tracer.Mainly from tungsten.
int generateLightPath(const Distribution1D *ptr, Sampler &sampler, int maxDepth, PathVertex *path);

int generateCameraPath();

int
randomWalk(const Scene &scene, Ray &ray, Sampler &sampler, int maxDepth, PathVertex *pVertex, bool adjoint, Float pdf);

Spectrum
connectPath(const Scene &scene, const Camera *camera, const PathVertex *lightPath, int ln, const PathVertex *cameraPath,
            int cn, float *misWeight, ivec2 *pRaster);

class BDPTIntegrator : public SamplerIntegrator {
public:
    BDPTIntegrator(std::shared_ptr<Camera> camera, std::shared_ptr<Sampler> sampler,
                   const Json &json) : SamplerIntegrator(camera, sampler, json){}

    void process(const Scene &scene, Sampler &sampler) override;

    vec3 integrate(const Ray &ray, const Scene &scene, Sampler &sampler) const override;

    void render(const Scene &scene) override;

protected:
    static inline int pyramidCount(int pathLength) {
        return ((pathLength + 1) * (pathLength + 2)) / 2 - 1;
    }

    static inline int pyramidIndex(int s, int t) {
        return pyramidCount(s + t - 2) + t - 1;
    }



    std::unique_ptr<Distribution1D> lightDistrib{nullptr};
    std::unique_ptr<std::map<const Light *,size_t>> lightIdx{};
    ImagePrid *imagePramid;
    std::vector<BdptTracer> _tracers;

    void saveOutPuts(const std::string &fileName,int spp);
};