#include "VolPathIntegrator.h"
#include "Bsdfs/Reflection.hpp"
#include "Mediums/Medium.hpp"


vec3 VolPathIntegrator::integrate(const Ray & ray, const Scene & scene, Sampler & sampler) const {
    Spectrum L(0);
    const Medium * medium;
    medium = _camera->_medium.get();
    VolumeEvent voulumeEvent;
    SurfaceEvent surfaceEvent;
    bool specularBounce = true;
    int bounce;
    Spectrum throughPut(1);
    Ray _ray(ray);
    for ( bounce = 0 ; bounce < maxBounces ; bounce ++ ) {
        std::optional < Intersection > its = scene.intersect(_ray);
        bool foundIntersection = its.has_value();
        //if(foundIntersection) return (its->Ng +1.f)/2.f; else return L;
        if ( ! foundIntersection && ! medium ) break;
        bool hitSurface = true;
        if ( medium ) {
//            return medium->sampleDistance(_ray, sampler, voulumeEvent);
            throughPut *= medium->sampleDistance(_ray, sampler, voulumeEvent);
            if ( isBlack(throughPut) ) break;
            hitSurface = voulumeEvent.exited;
        }
        if ( hitSurface ) {
//            return Spectrum(0);
            if ( specularBounce ) {
                if ( foundIntersection )
                    L += throughPut * its->Le(- _ray.d);
                else
                    for ( auto light: scene.lights ) {
                        if ( light->flags == int(LightFlags::Infinite) ) {
                            L += throughPut * light->Le(_ray);
                        }
                    }
            }
            surfaceEvent = makeLocalScatterEvent(& its.value());
            if ( its->bsdf->Pure(BSDF_FORWARD) ) {
                _ray = surfaceEvent.sctterRay(_ray.d);
            } else {
                if ( its->bsdf->MatchesFlags(BXDFType(BSDF_NO_SPECULAR)) && bounce < maxBounces - 1 ) {

                    L += throughPut * uniformSampleAllLights
                            (surfaceEvent, scene, sampler, medium);  //direct lighting
                }
                surfaceEvent.requestType = BSDF_ALL;
                Spectrum f = its->bsdf->sampleF(surfaceEvent, sampler.getNext2D(), false);
                if ( isBlack(f) || surfaceEvent.pdf == 0 )
                    break;
                BXDFType flags = surfaceEvent.sampleType;
                specularBounce = ( flags & BSDF_SPECULAR ) != 0;
                throughPut *= f / surfaceEvent.pdf;

                _ray = surfaceEvent.sctterRay(surfaceEvent.toWorld(surfaceEvent.wi));
            }
            medium = its->primitive->selectMedium(medium, dot(_ray.d, its->Ng) > 0);
        } else {
            //   return volumeUniformSampleOneLight(voulumeEvent,medium,scene,sampler,lightDistribution.get());
            L += throughPut *
                 volumeUniformSampleOneLight(voulumeEvent, medium, scene, sampler, lightDistribution.get());
            PhaseSample phaseSample;
            Spectrum p = voulumeEvent.phase->sampleP(_ray.d, sampler.getNext2D(), phaseSample);
            //   throughPut *= p / voulumeEvent.pdf;
            _ray = Ray(voulumeEvent.p, phaseSample.w);
        }
        //russian prob
        Float roulettePdf = std::max(throughPut.x, std::max(throughPut.y, throughPut.z));
        if ( bounce > 2 && roulettePdf < 0.1 ) {
            if ( sampler.getNext1D() < roulettePdf )
                throughPut /= roulettePdf;
            else
                break;
        }
    }
    return L;
}
