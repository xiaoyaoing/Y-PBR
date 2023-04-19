#include "scene.hpp"
#include "Ray/Ray.hpp"
#include "Camera/Camera.hpp"
#include "Sampler/Distrib.hpp"
#include "Sampler/Sampler.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"
class TraceHelper{
public:
    static bool sampleBSDF;
    static bool sampleLight;
    static Spectrum uniformSampleOneLight(SurfaceEvent & event, const Scene &scene,
                                   Sampler &sampler,
                                   const Distribution1D *lightDistrib = nullptr,
                                   const Medium * medium = nullptr);

    static Spectrum uniformSampleAllLights(SurfaceEvent & event, const Scene &scene,
                                    Sampler &sampler, const Medium * medium = nullptr);

    static Spectrum estimateDirect(SurfaceEvent & event, const vec2 & uShading, const Light & light, const vec2 & uLight,
                            const Scene & scene, Sampler & sampler, const Medium * medium = nullptr,
                            bool specular = false);

    static Spectrum volumeUniformSampleOneLight(VolumeEvent & event, const Medium * medium,const Scene &scene,
                                         Sampler &sampler,
                                         const Distribution1D *lightDistrib = nullptr
    ) ;

    static Spectrum volumeUniformSampleAllLights(VolumeEvent & event,const Medium * medium, const Scene &scene,
                                          Sampler &sampler);

    static Spectrum volumeEstimateDirect(VolumeEvent & event,const Medium * medium, const vec2 &uShading,
                                  const Light &light, const vec2 &uLight,
                                  const Scene &scene, Sampler &sampler
    ) ;


    static Spectrum evalLightDirect(const Scene & scene, const Light & light, Ray & ray,
                             const Medium * medium, Float * lightPdf);
    static Spectrum evalShadowDirect(const Scene & scene,Ray ray,const Medium * medium)  ;

    static std::unique_ptr<Distribution1D> computeLightPowerDistrib(const Scene & scene) ;

    static SurfaceEvent makeLocalScatterEvent(const Intersection * its);
};

bool russian(int depth, Sampler &sampler, Spectrum &beta);
