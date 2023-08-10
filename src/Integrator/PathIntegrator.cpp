#include "PathIntegrator.hpp"

#include "Bsdfs/Reflection.hpp"
#include "Common/Debug.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"
#include "Bsdfs/BSSRDF.hpp"
#include "TraceHelper.h"
#include <optional>

#include "iostream"
//2022/7/15
Spectrum PathIntegrator::integrate(const Ray &ray, const Scene &scene, Sampler &sampler) const {

    std::optional<Intersection> its;
    SurfaceEvent surfaceEvent;
    Spectrum beta(1.0);
    Spectrum L(0);

    int bounces = minBounces, maxDepth = maxBounces;
    bool specularBounce = true;
    Ray _ray(ray);

    std::string firstHitName;
    std::string Path;
    for (bounces = 0;; ++bounces) {

        its = scene.intersect(_ray);

        if (specularBounce && bounces>=minBounces) {
            if (its.has_value())
                L += beta * its->Le(-_ray.d);
            else
                for (auto light: scene.lights) {
                    if (light->flags == int(LightFlags::Infinite)) {
                        L += beta * light->Le(_ray);
                    }
                }
            if(hasNeg(L)){
                int k = 1;
            }
        }
        //  break;

        if (!its.has_value() || bounces >= maxDepth)
            break;
    //    return Spectrum(its->uv.x,its->uv.y,0);
        if (DebugConfig::OnlyShowNormal) {
            return (its->Ng + Spectrum(1.f)) / 2.f;
        }
        surfaceEvent = makeLocalScatterEvent(&its.value());
        its->bsdf->sampleF(surfaceEvent,sampler.getNext2D(),false);
        if (its->bsdf->Pure(BSDF_FORWARD)) {
            _ray = surfaceEvent.sctterRay(_ray.d);
        } else {
            if (!its->bsdf->Pure(BSDF_PURE_SPECULR) && bounces < maxDepth - 1) {
                Spectrum Ld = uniformSampleAllLights
                        (surfaceEvent, scene, sampler, nullptr);  //direct lighting
                if (DebugConfig::OnlyDirectLighting)
                    return Ld;
                 L += beta * Ld;
                if(hasNan(L)){
                    int k = 1;
                }

            }
            surfaceEvent.requestType = BSDF_ALL;
            Spectrum f = its->bsdf->sampleF(surfaceEvent, sampler.getNext2D(), false);
            if (isBlack(f) || surfaceEvent.pdf == 0)
                break;
            BXDFType flags = surfaceEvent.sampleType;
            specularBounce = (flags & BSDF_SPECULAR) != 0;
            beta *= f / surfaceEvent.pdf;

            if(hasNan(beta)){
                int k = 1;
            }
            _ray = surfaceEvent.sctterRay();


            if (its->bssrdf && (flags & BSDF_TRANSMISSION) ) {
                Intersection pi;
                Float pdf = 0;
                Spectrum s = its->bssrdf->sampleS(scene, sampler.getNext1D(), sampler.getNext2D(), surfaceEvent, &pi,
                                                  &pdf);
                if (isBlack(s) || pdf == 0)
                    break;
                surfaceEvent = makeLocalScatterEvent(&pi);
                beta *= s / pdf;
                L += beta * uniformSampleOneLight(surfaceEvent, scene, sampler, lightDistribution.get());
                surfaceEvent.requestType = BSDF_ALL;
                f = pi.bsdf->sampleF(surfaceEvent, sampler.getNext2D(), false);
                beta *= f / surfaceEvent.pdf;
                specularBounce = isSpecualr(surfaceEvent.sampleType);
                _ray = surfaceEvent.sctterRay();
            }
            if(russian(bounces,sampler,beta))
                break;
        }
    }
    if(hasNeg(L)){
        int k = 1;
    }
    //return vec3(bounces/10.f);
    return L;
}


void PathIntegrator::process(const Scene &scene, Sampler &sampler) {
    lightDistribution = CreateLightSampleDistribution(lightSampleStrategy, scene);
}
