#include "Integrator.hpp"
#include "iostream"
#include "Bsdfs/Reflection.hpp"

//static bool useMis = false;
static bool sampleBSDF = true;
static bool sampleLgiht = true;



Integrator::Integrator(nlohmann::json j) {

}


Spectrum
Integrator::EstimateDirect(SurfaceScatterEvent & event, const vec2 & uShading,
                           const Light & light, const vec2 & uLight,
                           const Scene & scene, Sampler & sampler, bool specular) const {

    BXDFType bsdfFlags =
            specular ? BSDF_ALL : BXDFType(BSDF_ALL & ~ BSDF_SPECULAR);
    event.sampleType = bsdfFlags;

    vec3 wi;
    Float lightPdf = 0, scatteringPdf = 0;
    VisibilityTester visibility;

    Spectrum Ld(0);
    //return Spectrum(17,12,4)/20.f;
    //sample light
    if(sampleLgiht)
    {
        Spectrum Li = light.Sample_Li(* ( event.its ), uLight, & wi, & lightPdf, & visibility);
        if ( ! isBlack(Li) && lightPdf != 0 ) {
            if ( visibility.Unoccluded(scene) ) {
              //  return Li;
                event.wi =  event.toLocal(wi);
                scatteringPdf = event.its->bsdf->Pdf(event);
                Spectrum f = event.its->bsdf->f(event) * abs(event.wi.z) ;
                if ( ! isBlack(f) ) {
                    if ( light.isDeltaLight() ) Ld += f * Li ;/// lightPdf;
                    else {
                        Float weight =
                                PowerHeuristic(lightPdf, scatteringPdf);
                        if ( !sampleBSDF ) weight = 1;
                        Ld += f * weight * Li /lightPdf;/// lightPdf;
                    }
                }
            }
        }
    }

    if(sampleBSDF)
    //sample bsdf
    {
        if ( ! light.isDeltaLight() ) {
            Spectrum f = event.its->bsdf->sampleF(event, uShading);
            scatteringPdf = event.pdf;
            f *= abs(event.wi.z);
         //   return event.toWorld(event.wi);
         //  return f;
            Float weight = 1;
            vec3 worldShadowRayDir = event.toWorld(event.wi);
          //  return (worldShadowRayDir+vec3(1))/2.f;
            Ray shaowRay(event.its->p + worldShadowRayDir  * Constant::EPSILON, worldShadowRayDir);
            //shaowRay.farT =100;
            std::optional < Intersection > its = scene.intersect(shaowRay);
            bool isLightInfinite = (light.flags & (int)LightFlags::Infinite);
            if ( (its.has_value() && its->primitive->areaLight.get() == & light) || isLightInfinite  ) {
                if(isLightInfinite){
                    if(its.has_value())
                    {
                          return Ld;
                    }
                    Intersection infinteIts;
                    its = std::make_optional(infinteIts);
                }
                its->w = shaowRay.d;
                if (!isSpecualr(event.sampleType) ) {

                    lightPdf = light.directPdf(its.value(), event.its->p);

//                    if ( lightPdf == 0 )
//                        return Ld;
                    weight = PowerHeuristic(scatteringPdf, lightPdf);
                    if(!sampleLgiht) weight = 1;
                }
            } else {
                return Ld;
            }
           // its->w = normalize(vec3(0.1));
            Spectrum Li = light.directLighting(its.value());
            Ld += f * weight * Li / scatteringPdf;
        }

    }
    return Ld;
}

Spectrum
Integrator::UniformSampleOneLight(SurfaceScatterEvent & event,
                                  const Scene & scene,
                                  Sampler & sampler,
                                  const Distribution1D * lightDistrib,
                                  bool handleMedia) const {

    //todo multiple lightPdf
    Float lightPdf;
    std::shared_ptr < Light > light;
    int lightNum;
    if ( lightDistrib ) {
        lightNum = lightDistrib->SampleDiscrete(sampler.getNext1D(), & lightPdf);
        if ( lightPdf == 0 ) return Spectrum(0.f);
    } else {
        int nLights = scene.lights.size();
        lightNum = std::min((int) ( sampler.getNext1D() * nLights ), nLights - 1);
        lightPdf = Float(1) / nLights;
    }
    light = scene.lights.at(lightNum);

    vec2 uScattering = sampler.getNext2D();
    vec2 uLight = sampler.getNext2D();
    return EstimateDirect(event, uScattering, * light, uLight,
                          scene, sampler, handleMedia) / lightPdf;
}
