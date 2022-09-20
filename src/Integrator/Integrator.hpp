#pragma  once

#include "../scene.hpp"
#include "../Common/math.hpp"
#include "../Sampler/Sampler.hpp"
#include "../Ray/Ray.hpp"
#include "../Sampler/Distrib.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"

#include <nlohmann/json.hpp>

struct LightSample{

};

class Integrator{
public:
   Integrator(nlohmann::json j);

   virtual vec3  integrate(const Ray& ray, const Scene& scene, Sampler& sampler) const =0;


   virtual void  Preprocess(const Scene &scene, Sampler & sampler )  =0;

    Spectrum UniformSampleOneLight(SurfaceScatterEvent & event, const Scene &scene,
                                   Sampler &sampler,
                                   const Distribution1D *lightDistrib,
                                   bool handleMedia = false) const;

    Spectrum UniformSampleAllLights(SurfaceScatterEvent & it, const Scene &scene, Sampler &sampler,
                                    const std::vector<int> &nLightSamples,
                                    bool handleMedia = false);

    Spectrum EstimateDirect(SurfaceScatterEvent & event, const vec2 &uShading,
                            const Light &light, const vec2 &uLight,
                            const Scene &scene, Sampler &sampler,
                            bool specular = false) const;

};