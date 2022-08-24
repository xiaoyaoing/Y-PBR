#include "PathIntegrator.hpp"
#include "optional"
#include "spdlog/spdlog.h"
//2022/7/15
Spectrum PathIntegrator::integrate(const Ray & ray, const Scene & scene, Sampler & sampler) const {
//    Intersection intersection =


//    vec3 lightPos(-20, 40, 20);
//    vec3 lightPower(2e3, 2.e3, 2e3);
//    vec3 dir = lightPos - intersection.pos;
//    return lightPower / length2(dir);

    std::optional <Intersection> its;
    Spectrum throughPut(1.0);
    Spectrum L(0);
    int bounces = 0, maxDepth = 10;
    Ray _ray(ray);
    Spectrum  tempL(0);
    for ( bounces = 0 ;; ++ bounces ) {
        if(throughPut.x<0 || throughPut.y<0 || throughPut.z<0){
            spdlog::info("invalid throughPut");
        }
        its = scene.intersect(_ray);

        if ( ! its.has_value() || bounces >= maxDepth ) break;

//        vec3 n = its->getNormal();
//        n += 1.0;
//        n/=2;
      //  return  n;

        L += throughPut * its->Le(- ray.d);   //hit emssive

        ///sample direct lighting for no specular bxdf
        its->wo = normalize(its->toLocal(- ray.d));
        if ( its->bsdf->NumComponents(BXDFType(BSDF_ALL & ~ BSDF_SPECULAR)) >
             0 ) {

            const Distribution1D *distrib = lightDistribution->Lookup(its->p);
            auto Ld = UniformSampleOneLight
                    (its.value(), scene, sampler,distrib);  //direct lighting
            if(bounces==0)
             tempL=L;
            L += throughPut * Ld;
            if( isnan(L.x)){
                int k=1;
            }
        }
    //   return L;
        Float pdf;
        BXDFType flags;
        vec3 wi;
        Spectrum f = its->bsdf->sampleF(its->wo, & wi, sampler.getNext2D(), & pdf,
                                        BSDF_ALL, &flags);
        throughPut *= f;

        ///Russian roulette to avoid Long distance path tracing(Unbiased estimation)
        Float roulettePdf = std::max(throughPut.x, std::max(throughPut.y, throughPut.z));
        if ( bounces > 4 && roulettePdf < 0.1 ) {
            if ( sampler.getNext1D() < roulettePdf )
                throughPut /= roulettePdf;
            else
                break;
        }

        _ray = Ray(its->p, its->toWorld(wi));
//         auto wo = its.shFrame.toLocal(-ray.d);
//
//         BsdfSample bsdfSample;
//         its->bsdf->sample(bsdfSample);
//
//         throughPut *= bsdfSample.val * cositem;
//         ray = Ray(its->pos,dir);
    }

    return L;
}

PathIntegrator::PathIntegrator(nlohmann::json j) : Integrator(j) ,
                                                   //lightSampleStrategy(j["lightSampleStrategy"])
                                                   lightSampleStrategy("uniform")
                                                   {
}

void PathIntegrator::Preprocess(const Scene & scene, Sampler & sampler) {
    lightDistribution = CreateLightSampleDistribution(lightSampleStrategy,scene);
}
