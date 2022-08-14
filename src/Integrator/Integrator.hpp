#pragma  once


#include <nlohmann/json.hpp>
#include "../scene.hpp"
#include "../Common/math.hpp"
#include "../Sampler/Sampler.hpp"
#include "../Ray/Ray.hpp"
#include "../Sampler/Distrib.hpp"
struct LightSample{

};

class Integrator{
public:
   Integrator(nlohmann::json j);

   virtual vec3  integrate(const Ray& ray, const Scene& scene, Sampler& sampler) const =0;

   virtual vec3  sampleDirectLight(const Scene & scene,const Intersection & intersection,Light sample);

   virtual void  Preprocess(const Scene &scene, Sampler & sampler )  =0;

    Spectrum UniformSampleOneLight(const Intersection &it, const Scene &scene,
                                   Sampler &sampler,
                                   const Distribution1D *lightDistrib,
                                   bool handleMedia = false) const;

    Spectrum UniformSampleAllLights(const Intersection &it, const Scene &scene,Sampler &sampler,
                                    const std::vector<int> &nLightSamples,
                                    bool handleMedia = false);

    Spectrum EstimateDirect(const Intersection &it, const vec2 &uShading,
                            const Light &light, const vec2 &uLight,
                            const Scene &scene, Sampler &sampler,
                            bool specular = false) const;

};