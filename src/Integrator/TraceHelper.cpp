#include "TraceHelper.h"
#include "Bsdfs/Reflection.hpp"
#include "Common/ProgressReporter.h"
#include "Sampler/Sampler.hpp"
#include "Common/Parallel.h"
#include "Camera/Camera.hpp"

#include <thread>
#include <iostream>
#include "Texture/BitMapTexture.hpp"
#include "Mediums/Medium.hpp"

bool TraceHelper::sampleBSDF  = true;
bool TraceHelper::sampleLight = true;

bool sampleLight = TraceHelper::sampleLight;
bool sampleBSDF  = TraceHelper::sampleBSDF;

Spectrum estimateDirect(SurfaceEvent& event, const vec2& uShading, const Light& light, const vec2& uLight, const Scene& scene, Sampler& sampler, const Medium* medium, bool specular) {
    const BSDF* bsdf = event.its->bsdf;

    event.requestType =
        specular ? BSDF_ALL : BSDF_NO_SPECULAR;

    //    if ( ! event.its->bsdf->MatchesFlags(event.requestType) )
    //        return Spectrum(0);
    //    assert(!event.its->bsdf->HasFlag(BSDF_TRANSMISSION));
    vec3             wi;
    Float            lightPdf = 0, scatteringPdf = 0;
    VisibilityTester visibility;
    Spectrum         Ld(0);
    if (sampleLight) {
        Float    distance;
        Spectrum Li = light.sampleLi(event.its->p, uLight, &wi, &lightPdf, &distance);
        Ray      ray(event.its->p, wi, Constant::EPSILON, distance - Constant::EPSILON);
        Li *= evalShadowDirect(scene, ray, medium);
        if (!isBlack(Li) && lightPdf != 0) {

            auto its      = scene.intersect(ray);
            event.wi      = event.toLocal(wi);
            scatteringPdf = bsdf->Pdf(event);
            Spectrum f    = bsdf->f(event, false);

            if (!isBlack(f)) {
                if (light.isDeltaLight())
                    Ld += f * Li;/// lightPdf;
                else {
                    Float weight =
                        PowerHeuristic(lightPdf, scatteringPdf);
                    if (!sampleBSDF) weight = 1;
                    Ld += f * weight * Li / lightPdf;
                    if (hasNan(Ld)) {
                        1;
                    }
                }
            }
        }
    }
    if (sampleBSDF) {
        if (!light.isDeltaLight()) {
            // return event.its->Ns;
            Spectrum f    = bsdf->sampleF(event, uShading, false);
            scatteringPdf = event.pdf;
            if (!isBlack(f) && scatteringPdf != 0) {
                vec3 worldShadowRayDir = event.toWorld(event.wi);
                Ray  shaowRay(event.its->p, worldShadowRayDir);
                Spectrum Li = evalLightDirect(scene, light, shaowRay, medium, &lightPdf);
                if (lightPdf == 0 || isBlack(Li)) return Ld;
                Float weight = 1;
                if (!isSpecualr(event.sampleType)) {
                    weight = PowerHeuristic(scatteringPdf, lightPdf);
                    if (!sampleLight) weight = 1;
                }

                Ld += f * weight * Li / scatteringPdf;
                if (hasNan(Ld)) {
                    bsdf->f(event, false);
                    1;
                }
            }
        }
    }
    if (hasNan(Ld)) {
        DebugBreak();
    }
    // assert(isBlack(Ld));
    return Ld;
}

Spectrum
uniformSampleOneLight(SurfaceEvent&         event,
                      const Scene&          scene,
                      Sampler&              sampler,
                      const Distribution1D* lightDistrib,
                      const Medium*         medium) {

    //todo multiple lightPdf
    Float                  lightPdf;
    std::shared_ptr<Light> light;
    int                    lightNum;
    if (lightDistrib) {
        lightNum = lightDistrib->SampleDiscrete(sampler.getNext1D(), &lightPdf);
        if (lightPdf == 0) return Spectrum(0.f);
    } else {
        int nLights = scene.lights.size();
        lightNum    = std::min((int)(sampler.getNext1D() * nLights), nLights - 1);
        lightPdf    = Float(1) / nLights;
    }
    light = scene.lights.at(lightNum);

    vec2 uScattering = sampler.getNext2D();
    vec2 uLight      = sampler.getNext2D();
    return estimateDirect(event, uScattering, *light, uLight, scene, sampler, medium, false) / lightPdf;
}

Spectrum uniformSampleAllLights(SurfaceEvent& event, const Scene& scene, Sampler& sampler, const Medium* medium) {
    Spectrum L(0);
    for (auto& light : scene.lights) {
        L += estimateDirect(event, sampler.getNext2D(), *light, sampler.getNext2D(), scene, sampler, medium, false);
    }
    return L;
}

SurfaceEvent makeLocalScatterEvent(const Intersection* its) {
    SurfaceEvent event;
    Frame        frame = its->primitive->setTangentFrame(its);

    bool enableTwoSideShading = true;

    bool hitBackSide    = dot(its->w, its->Ng) > 0;
    bool isTransmissive = its->bsdf->HasFlag(BSDF_TRANSMISSION);
    bool flippedFrame   = false;
    //todo add this to config class
    if (enableTwoSideShading && hitBackSide && !isTransmissive) {
        frame.n       = -frame.n;
        frame.tangent = -frame.tangent;
        flippedFrame  = true;
    }
    event.frame        = frame;
    event.its          = its;
    event.wo           = event.toLocal(-its->w);
    event.flippedFrame = flippedFrame;
    event.requestType  = BSDF_ALL;
    return event;
}

std::unique_ptr<Distribution1D> computeLightPowerDistrib(const Scene& scene) {
    int    nLights = scene.lights.size();
    Float* power   = new Float[nLights];
    for (int i = 0; i < nLights; i++) {
        power[i] = luminace(scene.lights[i]->Power());
    }
    return std::make_unique<Distribution1D>(power, nLights);
}

Spectrum
evalLightDirect(const Scene& scene, const Light& light, Ray& ray, const Medium* medium, Float* lightPdf) {
    if (light.isDeltaLight())
        return Spectrum(0);
    std::optional<Intersection> lightIts = light.intersect(ray);
    if (!lightIts) return Spectrum();
    Spectrum L;
    if (light.flags & int(LightFlags::Area)) L = lightIts->Le(-ray.d);
    if (light.flags & int(LightFlags::Infinite)) L = light.Le(ray);
    //avoid self shadow
    if (isBlack(L)) {
        if(lightPdf) *lightPdf = 0;
        return Spectrum(0);
    }
    if (lightPdf) *lightPdf = light.PdfLi(lightIts.value(), ray.o);
    return L * evalShadowDirect(scene, ray, medium);
}

Spectrum evalShadowDirect(const Scene& scene, Ray ray, const Medium* medium) {
    if (!medium) {
        return scene.intersectP(ray) ? Spectrum(0) : Spectrum(1);
    }
    Spectrum                    Tr(1);
    std::optional<Intersection> its;
    Float                       tHit = ray.farT;
    while (tHit > 0 && medium) {
        its         = scene.intersect(ray);
        bool didHit = its.has_value();
        if (didHit && !its->bsdf->Pure(BSDF_FORWARD))
            return Spectrum();
        Tr *= medium->TR(ray);
        if (!its)
            return Tr;
        medium = its->primitive->selectMedium(medium, dot(ray.d, its->Ng) > 0);
        tHit -= ray.farT;
        ray = Ray(its->p, ray.d, Constant::EPSILON, tHit);
    }
    return Tr;
}

Spectrum volumeUniformSampleOneLight(VolumeEvent& event, const Medium* medium, const Scene& scene, Sampler& sampler, const Distribution1D* lightDistrib) {
    Float                  lightPdf;
    std::shared_ptr<Light> light;
    int                    lightNum;
    if (lightDistrib) {
        lightNum = lightDistrib->SampleDiscrete(sampler.getNext1D(), &lightPdf);
        if (lightPdf == 0) return Spectrum(0.f);
    } else {
        int nLights = scene.lights.size();
        lightNum    = std::min((int)(sampler.getNext1D() * nLights), nLights - 1);
        lightPdf    = Float(1) / nLights;
    }
    lightNum = 0;
    lightPdf = 1;
    light    = scene.lights.at(lightNum);

    vec2 uScattering = sampler.getNext2D();
    vec2 uLight      = sampler.getNext2D();
    return volumeEstimateDirect(event, medium, uScattering, *light, uLight, scene, sampler) / lightPdf;
}

Spectrum volumeUniformSampleAllLights(VolumeEvent& event, const Medium* medium, const Scene& scene, Sampler& sampler) {
    return Spectrum();
}

Spectrum
volumeEstimateDirect(VolumeEvent& event, const Medium* medium, const vec2& uShading, const Light& light, const vec2& uLight, const Scene& scene, Sampler& sampler) {
    Spectrum Ld(0), Li;
    vec3     wi;
    Float    lightPdf, distance, scatteringPdf;
    ///light sample
    Li = light.sampleLi(event.p, uLight, &wi, &lightPdf, &distance);
    Ray ray(event.p, wi, Constant::EPSILON, distance - Constant::EPSILON);
    // return Spectrum(distance);
    Spectrum p = event.phase->p(event.rayDir, wi);
    if (!isBlack(p) && !isBlack((Li = Li * evalShadowDirect(scene, ray, medium)))) {
        Float weight  = 1;
        scatteringPdf = event.phase->pdf(-event.rayDir, wi);
        if (!light.isDeltaLight()) weight = PowerHeuristic(lightPdf, scatteringPdf);
        Ld += weight * p * Li / lightPdf;
        //   return Ld;
        // return vec3(lightPdf/100.f);
        //      return Li/lightPdf;
    }
    ///phase sample

    if (!light.isDeltaLight()) {
        PhaseSample phaseSample;
        Spectrum    f = event.phase->sampleP(event.rayDir, uShading, phaseSample);
        scatteringPdf = phaseSample.pdf;
        if (!isBlack(f) && scatteringPdf) {
            Ray shaowRay(event.p, phaseSample.w, Constant::EPSILON);
            Li = evalLightDirect(scene, light, shaowRay, medium, &lightPdf);
            if (lightPdf == 0) return Ld;
            Float weight = PowerHeuristic(scatteringPdf, lightPdf);
            if (!isBlack(Li)) Ld += f * weight * Li / scatteringPdf;
            if (luminace(Ld) > 10) {
                DebugBreak();
            }
        }
    }
    return Ld;
}

bool russian(int depth, Sampler& sampler, Spectrum beta) {
    float pdf = max(beta);
    if (depth < 2 || pdf > 0.1)
        return false;
    if (sampler.getNext1D() < pdf) {
        beta /= pdf;
        return false;
    }
    return true;
}