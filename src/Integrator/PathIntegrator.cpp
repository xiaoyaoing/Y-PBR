#include "PathIntegrator.hpp"

#include "Bsdfs/Reflection.hpp"
#include "Common/Debug.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"

#include <optional>
#include <spdlog/spdlog.h>

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
    bool specularBounce = true;
    Ray _ray(ray);

    std::string firstHitName ;
    std::string Path;
    for ( bounces = 0 ;; ++ bounces ) {
        if(abs(_ray.d.x- -0.027210)<0.01){
            int k=1;
        }
        its = scene.intersect(_ray);

        if(specularBounce)
        {
            if(its.has_value())
                L += throughPut * its->Le(-_ray.d);
            else for(auto light :scene.lights){
                if(light->flags == int(LightFlags::Infinite)){
                    L +=light->environmentLighting(_ray);
                }
            }
        }  //hit emssive
        if ( ! its.has_value() || bounces >= maxDepth ) break;

//        char raydirStr[128];
//        sprintf(raydirStr,"%1f %1f %1f",_ray.d.x,_ray.d.y,_ray.d.z);
//        if(bounces == 0)
//        firstHitName = its->bsdf->name;
//        Path+=" "+its->bsdf->name+" "+std::string(raydirStr);
        if(DebugConfig::OnlyShowNormal){
                return its->Ng;
        }



        SurfaceScatterEvent localScatter = makeLocalScatterEvent(&its.value());

        ///sample direct lighting for no specular bxdf
        if ( its->bsdf->MatchesFlags(BXDFType(BSDF_ALL & ~ BSDF_SPECULAR)) ) {

            const Distribution1D *distrib = lightDistribution->Lookup(its->p);
            auto Ld = UniformSampleOneLight
                    (localScatter, scene, sampler,distrib);  //direct lighting
            if(DebugConfig::OnlyDirectLighting)
                return Ld;
            L += throughPut * Ld;
        }

        if(DebugConfig::OnlyIndirectLighting && bounces==0){
            L = Spectrum (0.f);
        }


        //   return L;

        Spectrum  f =its->bsdf->sampleF(localScatter,sampler.getNext2D());
        if( isBlack(f) || localScatter.pdf==0)
            break;

        BXDFType flags = localScatter.sampleType;
        specularBounce = (flags & BSDF_SPECULAR)!=0;

        vec3 newDir = localScatter.toWorld(localScatter.wi);
      //  newDir = vec3(sampler.getNext1D(),sampler.getNext1D(),sampler.getNext1D());
        if(specularBounce) throughPut *= f;
        else throughPut *= f * absDot(newDir,its->Ns) / localScatter.pdf;

        if( isnan(throughPut.x)){
            int k =1;
        }


        if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)){
           Float  eta = its->bsdf->eta();
           Float  etaScale = (dot(newDir, its->Ng) > 0) ? (eta * eta) : 1 / (eta * eta);
           throughPut *= etaScale;
        }

        _ray = Ray(its->p, newDir,Constant::EPSILON);

       auto  tempIts = scene.intersect(_ray);
        ///Russian roulette to avoid Long distance path tracing(Unbiased estimation)
        Float roulettePdf = std::max(throughPut.x, std::max(throughPut.y, throughPut.z));
        if ( bounces > 4 && roulettePdf < 0.1 ) {
            if ( sampler.getNext1D() < roulettePdf )
                throughPut /= roulettePdf;
            else
                break;
        }

    }

    if( isnan(L.x)){
        int k=1;
    }
//    if()
   // spdlog::info("{0}", toColorStr(L));

    if( firstHitName=="water")  {
      //  spdlog::info("bounce count {0} {1} {2}  ",bounces, toColorStr(L),Path);
        //L*=10;
    }
    return L;
}

PathIntegrator::PathIntegrator(nlohmann::json j) : AbstractPathTracer(j) ,
                                                   //lightSampleStrategy(j["lightSampleStrategy"])
                                                   lightSampleStrategy("uniform")
                                                   {
}

void PathIntegrator::Preprocess(const Scene & scene, Sampler & sampler) {
    lightDistribution = CreateLightSampleDistribution(lightSampleStrategy,scene);
}
