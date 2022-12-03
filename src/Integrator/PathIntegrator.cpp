#include "PathIntegrator.hpp"

#include "Bsdfs/Reflection.hpp"
#include "Common/Debug.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"

#include <optional>
#include <spdlog/spdlog.h>

//2022/7/15
Spectrum PathIntegrator::integrate(const Ray & ray, const Scene & scene, Sampler & sampler) const {
    std::optional < Intersection > its;
    SurfaceEvent surfaceEvent;
    Spectrum throughput(1.0);
    Spectrum L(0);

    int bounces = 0, maxDepth = maxBounces;
    bool specularBounce = true;
    Ray _ray(ray);

    std::string firstHitName;
    std::string Path;
    for ( bounces = 0 ;; ++ bounces ) {

        its = scene.intersect(_ray);
        if ( specularBounce) {
            if ( its.has_value() )
                L += throughput * its->Le(- _ray.d);
            else
                for ( auto light: scene.lights ) {
                    if ( light->flags == int(LightFlags::Infinite) ) {
                         L += throughput * light->Le(_ray);
                    }
                }
        }
      //  break;

        if ( ! its.has_value() || bounces >= maxDepth )
            break;
        if ( DebugConfig::OnlyShowNormal ) {
            return ( its->Ng + Spectrum(1.f) ) / 2.f;
        }
        surfaceEvent = makeLocalScatterEvent(& its.value());

        if ( its->bsdf->Pure(BSDF_FORWARD) ) {
            _ray = surfaceEvent.sctterRay(_ray.d);
        } else
            ///sample direct lighting for no specular bxdf
        {
            if ( its->bsdf->MatchesFlags(BXDFType(BSDF_NO_SPECULAR)) && bounces < maxDepth - 1 ) {

                Spectrum Ld = uniformSampleAllLights
                        (surfaceEvent, scene, sampler, nullptr);  //direct lighting
                if ( DebugConfig::OnlyDirectLighting )
                    return Ld;
                L += throughput * Ld;
                if ( hasNan(L) ) {
                    int k = 1;
                }
            }

            {
                if ( DebugConfig::OnlyDirectLighting )
                    return L;
                if ( DebugConfig::OnlyIndirectLighting && bounces == 0 ) {
                    L = Spectrum(0.f);
                }
            }

            surfaceEvent.requestType = BSDF_ALL;
            Spectrum f = its->bsdf->sampleF(surfaceEvent, sampler.getNext2D(), false);
            if ( isBlack(f) || surfaceEvent.pdf == 0 )
                break;

            BXDFType flags = surfaceEvent.sampleType;
            specularBounce = ( flags & BSDF_SPECULAR ) != 0;

            throughput *= f / surfaceEvent.pdf;
//            if(bounces == 1)
//            return (surfaceEvent.toWorld(surfaceEvent.wi)+vec3(1.f))/2.f;

            Float roulettePdf = std::max(throughput.x, std::max(throughput.y, throughput.z));
            if ( bounces > 2 && roulettePdf < 0.1 ) {
                if ( sampler.getNext1D() < roulettePdf )
                    throughput /= roulettePdf;
                else
                    break;
            }
            _ray = surfaceEvent.sctterRay(surfaceEvent.toWorld(surfaceEvent.wi));
        }
    }
    //return throughput;
    return L;
}


void PathIntegrator::process(const Scene & scene, Sampler & sampler) {
    lightDistribution = CreateLightSampleDistribution(lightSampleStrategy, scene);
}
