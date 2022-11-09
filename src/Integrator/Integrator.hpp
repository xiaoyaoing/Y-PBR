#pragma  once

#include "scene.hpp"
#include "Ray/Ray.hpp"
#include "Camera/Camera.hpp"
#include "Sampler/Distrib.hpp"
#include "Sampler/Sampler.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"




//class Sampler;

//class Image;


class Integrator{
public:
  // Integrator(Json j);
   virtual void render(const Scene & scene) const = 0 ;

   virtual void  process(const Scene &scene, Sampler & sampler )  =0;

    Spectrum uniformSampleOneLight(SurfaceScatterEvent & event, const Scene &scene,
                                   Sampler &sampler,
                                   const Distribution1D *lightDistrib = nullptr,
                                   bool handleMedia = false) const;

    Spectrum uniformSampleAllLights(SurfaceScatterEvent & event, const Scene &scene,
                                    Sampler &sampler,bool handleMedia = false) const ;

    Spectrum estimateDirect(SurfaceScatterEvent & event, const vec2 &uShading,
                            const Light &light, const vec2 &uLight,
                            const Scene &scene, Sampler &sampler,
                            bool specular = false) const;

    SurfaceScatterEvent makeLocalScatterEvent(const Intersection * its) const ;
    std::unique_ptr<Distribution1D> computeLightPowerDistrib(const Scene & scene) const;
};

//Integrator Based Sampler
class SamplerIntegrator : public  Integrator{

public:
    SamplerIntegrator(std::shared_ptr<Camera> camera,std::shared_ptr<Sampler> sampler):_camera(camera),_sampler(sampler){}
    void render(const Scene & scene) const override;
    virtual vec3  integrate(const Ray& ray, const Scene& scene, Sampler& sampler) const = 0 ;
private:
    std::shared_ptr<Camera> _camera;
    std::shared_ptr<Sampler> _sampler;
};