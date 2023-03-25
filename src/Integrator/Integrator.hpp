#pragma  once

#include "scene.hpp"
#include "Ray/Ray.hpp"
#include "Camera/Camera.hpp"
#include "Sampler/Distrib.hpp"
#include "Sampler/Sampler.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"

class Medium;

//class Sampler;

//class Image;


class Integrator{
public:
   Integrator(const Json & json){
       minBounces = getOptional(json,"min_bounces",0);
       maxBounces = getOptional(json,"max_bounces",8);
   }
   virtual void render(const Scene & scene) const = 0 ;

   virtual void  process(const Scene &scene, Sampler & sampler )  =0;

    Spectrum uniformSampleOneLight(SurfaceEvent & event, const Scene &scene,
                                   Sampler &sampler,
                                   const Distribution1D *lightDistrib = nullptr,
                                   const Medium * medium = nullptr) const;

    Spectrum uniformSampleAllLights(SurfaceEvent & event, const Scene &scene,
                                    Sampler &sampler, const Medium * medium = nullptr) const ;

    Spectrum estimateDirect(SurfaceEvent & event, const vec2 & uShading, const Light & light, const vec2 & uLight,
                            const Scene & scene, Sampler & sampler, const Medium * medium = nullptr,
                            bool specular = false) const;

    Spectrum volumeUniformSampleOneLight(VolumeEvent & event, const Medium * medium,const Scene &scene,
                                   Sampler &sampler,
                                   const Distribution1D *lightDistrib = nullptr
                                 ) const;

    Spectrum volumeUniformSampleAllLights(VolumeEvent & event,const Medium * medium, const Scene &scene,
                                    Sampler &sampler) const ;

    Spectrum volumeEstimateDirect(VolumeEvent & event,const Medium * medium, const vec2 &uShading,
                            const Light &light, const vec2 &uLight,
                            const Scene &scene, Sampler &sampler
                            ) const;


    Spectrum evalLightDirect(const Scene & scene, const Light & light, Ray & ray,
                             const Medium * medium, Float * lightPdf) const;
    Spectrum evalShadowDirect(const Scene & scene,Ray ray,const Medium * medium) const ;

    std::unique_ptr<Distribution1D> computeLightPowerDistrib(const Scene & scene) const;
protected :
    int minBounces,maxBounces;
};

//Integrator Based Sampler
class SamplerIntegrator : public  Integrator{

public:
    SamplerIntegrator(std::shared_ptr<Camera> camera,std::shared_ptr<Sampler> sampler,const Json & json):
    _camera(camera),_sampler(sampler), Integrator(json){}
    void render(const Scene & scene) const override;
    virtual vec3  integrate(const Ray& ray, const Scene& scene, Sampler& sampler) const = 0 ;

public:
    std::shared_ptr<Camera> _camera;
    std::shared_ptr<Sampler> _sampler;
};

SurfaceEvent makeLocalScatterEvent(const Intersection * its);
bool russian(int depth, Sampler &sampler, vec3 &beta);