#include "PathIntegrator.hpp"
#include "optional"
#include "spdlog/spdlog.h"
#include "Common/Debug.hpp"
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
       if(firstHitName=="IceAir" && bounces==1){
           int k=1;
       }
        its = scene.intersect(_ray);

        if ( ! its.has_value() || bounces >= maxDepth ) break;


        Path +=its->bsdf->name+" ";
        if(Path =="IceAir IceAir "){
            int k=1;
        }

        if(DebugConfig::OnlyShowNormal){
           return (its->getNormal()+1.0f)/2.f;
        }

        if(bounces==0){
            firstHitName = its->bsdf->name;
            if(firstHitName=="IceAir"){
                int k=1;
            }
        }

        if(its->bsdf->name == "lighting" && bounces!=0){
            int k=1;
        }

        if(specularBounce)
        L += throughPut * its->Le(- ray.d);   //hit emssive

        if(!its->bsdf->NumComponents(BSDF_TRANSMISSION) && dot(its->getNormal(),_ray.d )>0 )
        {
            its->shFrame. n = - its->shFrame.n;
//            its->setNormal(-its->getNormal());
            its->shFrame. s = - its->shFrame.s;
        }

        ///sample direct lighting for no specular bxdf
        its->wo = normalize(its->toLocal(- _ray.d));
        if ( its->bsdf->NumComponents(BXDFType(BSDF_ALL & ~ BSDF_SPECULAR)) >
             0 ) {

            const Distribution1D *distrib = lightDistribution->Lookup(its->p);
            auto Ld = UniformSampleOneLight
                    (its.value(), scene, sampler,distrib);  //direct lighting

            L += throughPut * Ld;

            if(DebugConfig::OnlyDirectLighting){
                return L;
            }
        }

        if(DebugConfig::OnlyIndirectLighting && bounces==0){
            L = Spectrum (0.f);
        }
    //   return L;
        Float pdf;
        BXDFType flags;
        vec3 wi;
        Spectrum f = its->bsdf->sampleF(its->wo, & wi, sampler.getNext2D(), & pdf,BSDF_ALL, &flags);

        specularBounce = (flags & BSDF_SPECULAR)!=0;

        vec3 newDir = its->toWorld(wi);
        if(specularBounce) throughPut *= f;
        else throughPut *= f * absDot(newDir,its->getNormal()) / pdf;

        if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)){
           Float  eta = 1.31;
           Float  etaScale = (dot(newDir, its->getNormal()) > 0) ? (eta * eta) : 1 / (eta * eta);
           throughPut *= etaScale;
        }

       //throughPut*=f;
       // spdlog::info(toColorStr(wi));

        _ray = Ray(its->p+newDir * Constant::EPSILON, newDir);

       auto  tempIts = scene.intersect(_ray);
        if(tempIts.has_value() && tempIts->primitive->bsdf->name=="IceAir" && bounces==0){
            int k=1;
        }

        ///Russian roulette to avoid Long distance path tracing(Unbiased estimation)
        Float roulettePdf = std::max(throughPut.x, std::max(throughPut.y, throughPut.z));
        if ( bounces > 4 && roulettePdf < 0.1 ) {
            if ( sampler.getNext1D() < roulettePdf )
                throughPut /= roulettePdf;
            else
                break;
        }

    }
//    if()
   // spdlog::info("{0}", toColorStr(L));

//    if( firstHitName=="IceAir")  {
//        spdlog::info("bounce count {0} {1} {2}  ",bounces, toColorStr(L),Path);
//    }
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
