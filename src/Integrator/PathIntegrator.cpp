#include "PathIntegrator.hpp"

#include "Bsdfs/Reflection.hpp"
#include "Common/Debug.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"

#include <optional>
#include <spdlog/spdlog.h>

//2022/7/15
Spectrum PathIntegrator::integrate(const Ray & ray, const Scene & scene, Sampler & sampler) const {

    std::optional < Intersection > its;
    Spectrum throughPut(1.0);
    Spectrum L(0);
    int bounces = 0, maxDepth = 10;
    bool specularBounce = true;
    Ray _ray(ray);

    std::string firstHitName;
    std::string Path;
    for ( bounces = 0 ;; ++ bounces ) {
        its = scene.intersect(_ray);
        if ( specularBounce ) {
            if ( its.has_value() )
                L += throughPut * its->Le(- _ray.d);
            else
                for ( auto light: scene.lights ) {
                    if ( light->flags == int(LightFlags::Infinite) ) {
                        L += throughPut * light->environmentLighting(_ray);
                    }
                }

            if( hasNan(L)){
                int k=1;
            }
        }  //hit emssive
        if ( ! its.has_value() || bounces >= maxDepth ) break;

        if ( DebugConfig::OnlyShowNormal ) {
            return its->Ng;
        }

        SurfaceScatterEvent localScatter = makeLocalScatterEvent(& its.value());

        ///sample direct lighting for no specular bxdf
        if ( its->bsdf->MatchesFlags(BXDFType(BSDF_ALL & ~ BSDF_SPECULAR))) {

            const Distribution1D * distrib = lightDistribution->Lookup(its->p);
            Spectrum Ld = uniformSampleOneLight
                    (localScatter, scene, sampler, distrib);  //direct lighting
            if ( DebugConfig::OnlyDirectLighting )
                return Ld;
            L += throughPut * Ld;
            if( hasNan(L)){
                int k=1;
            }
        }

        if ( DebugConfig::OnlyIndirectLighting && bounces == 0 ) {
            L = Spectrum(0.f);
        }

        Spectrum f = its->bsdf->sampleF(localScatter, sampler.getNext2D());
        if ( isBlack(f) || localScatter.pdf == 0 )
            break;

        BXDFType flags = localScatter.sampleType;
        specularBounce = ( flags & BSDF_SPECULAR ) != 0;

        vec3 newDir = localScatter.toWorld(localScatter.wi);

        throughPut *= f * absDot(newDir, its->Ns) / localScatter.pdf;


        if ( ( flags & BSDF_SPECULAR ) && ( flags & BSDF_TRANSMISSION ) ) {
            Float eta = its->bsdf->eta();
            Float etaScale = ( dot(newDir, its->Ng) > 0 ) ? ( eta * eta ) : 1 / ( eta * eta );
            throughPut *= etaScale;
        }
        else if(bounces>1){
            int k = 1;
        }

        _ray = Ray(its->p, newDir, Constant::EPSILON);

        std::optional<Intersection> tempIts = scene.intersect(_ray);
        ///Russian roulette to avoid Long distance path tracing(Unbiased estimation)
        Float roulettePdf = std::max(throughPut.x, std::max(throughPut.y, throughPut.z));
        if ( bounces > 128 && roulettePdf < 0.1 ) {
            if ( sampler.getNext1D() < roulettePdf )
                throughPut /= roulettePdf;
            else
                break;
        }
    }
    return L;
}

//PathIntegrator::PathIntegrator(Json j) : AbstractPathTracer(j),
//        //lightSampleStrategy(j["lightSampleStrategy"])
//                                                   lightSampleStrategy("uniform") {
//}

void PathIntegrator::process(const Scene & scene, Sampler & sampler) {
    lightDistribution = CreateLightSampleDistribution(lightSampleStrategy, scene);
}
