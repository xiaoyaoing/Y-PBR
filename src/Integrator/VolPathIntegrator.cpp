#include "VolPathIntegrator.h"
#include "Bsdfs/Reflection.hpp"
#include "Mediums/Medium.hpp"

void VolPathIntegrator::process(const Scene & scene, Sampler & sampler) {

}

vec3 VolPathIntegrator::integrate(const Ray & ray, const Scene & scene, Sampler & sampler) const {
    Spectrum L(0);
    const Medium * medium;
    medium = _camera->_medium.get();
    VolumeEvent voulumeEvent;
    SurfaceEvent surfaceEvent;
    bool specularBounce = true;
    int bounce;
    int maxBounce = 8;
    Spectrum throughPut(1);
    Ray _ray(ray);
    for ( bounce = 0 ; bounce < maxBounce ; bounce ++ ) {
        std::optional < Intersection > its = scene.intersect(_ray);
        bool foundIntersection = its.has_value();
        if ( ! foundIntersection && ! medium ) break;
        bool hitSurface = true;
        if ( medium ) {
            throughPut *= medium->sampleDistance(ray, sampler, voulumeEvent);
            if ( isBlack(throughPut) ) break;
            hitSurface = voulumeEvent.exited;
        }
        if ( hitSurface ) {
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
                if ( its->bsdf->MatchesFlags(BXDFType(BSDF_NO_SPECULAR)) && bounce < maxBounce - 1 ) {

                    Spectrum Ld = uniformSampleAllLights
                            (surfaceEvent, scene, sampler, lightDistribution.get());  //direct lighting
                    L += throughPut * Ld;
                }
                surfaceEvent.requestType = BSDF_ALL;
                Spectrum f = its->bsdf->sampleF(surfaceEvent, sampler.getNext2D(), false);
                if ( isBlack(f) || surfaceEvent.pdf == 0 )
                    break;

                BXDFType flags = surfaceEvent.sampleType;
                specularBounce = ( flags & BSDF_SPECULAR ) != 0;

                throughPut *= f / surfaceEvent.pdf;

                Float roulettePdf = std::max(throughPut.x, std::max(throughPut.y, throughPut.z));
                if ( bounce > 2 && roulettePdf < 0.1 ) {
                    if ( sampler.getNext1D() < roulettePdf )
                        throughPut /= roulettePdf;
                    else
                        break;
                }
                _ray = surfaceEvent.sctterRay(surfaceEvent.toWorld(surfaceEvent.wi));
            }
            medium = its->primitive->selectMedium(medium,dot(-_ray.d,its->Ng)>0);
        } else {
           // L +=
        }
    }

    return vec3();
}
