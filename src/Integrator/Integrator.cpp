#include "Integrator.hpp"
#include "Bsdfs/Reflection.hpp"
#include "Common/ProgressReporter.h"
#include "Sampler/Sampler.hpp"
#include "Common/Parallel.h"
#include "Camera/Camera.hpp"
#include "PathIntegrator.hpp"
#include <thread>
#include <iostream>

#include "Texture/BitMapTexture.hpp"

#include "Mediums/Medium.hpp"

static bool sampleBSDF = true;
static bool sampleLgiht = true;


//Integrator::Integrator(Json j) {
//
//}


Spectrum
Integrator::estimateDirect(SurfaceEvent & event, const vec2 & uShading, const Light & light, const vec2 & uLight,
                           const Scene & scene, Sampler & sampler, const Medium * medium, bool specular) const {
    const BSDF * bsdf = event.its->bsdf;

    //pure specular case
    if ( ! bsdf->MatchesFlags(BSDF_NO_SPECULAR) )
        return Spectrum();


    event.requestType =
            specular ? BSDF_ALL : BSDF_NO_SPECULAR;
    vec3 wi;
    Float lightPdf = 0, scatteringPdf = 0;
    VisibilityTester visibility;

    Spectrum Ld(0);
    if ( sampleLgiht ) {
        Float distance;
        Spectrum Li = light.sampleLi(event.its->p, uLight, & wi, & lightPdf, & distance);
        // return Spectrum(lightPdf);

        Ray ray(event.its->p, wi, Constant::EPSILON, distance - Constant::EPSILON);
        Li *= evalShadowDirect(scene, ray, medium);
        if ( ! isBlack(Li) && lightPdf != 0 ) {
            event.wi = event.toLocal(wi);
            scatteringPdf = bsdf->Pdf(event);
            Spectrum f = bsdf->f(event, false);
            if ( ! isBlack(f) ) {
                if ( light.isDeltaLight() ) Ld += f * Li;/// lightPdf;
                else {
                    Float weight =
                            PowerHeuristic(lightPdf, scatteringPdf);
                    if ( ! sampleBSDF ) weight = 1;
                    Ld += f * weight * Li / lightPdf;
                    if ( hasNan(Ld) ) {
                        bsdf->f(event, false);
                        int k = 1;
                    }
                }
            }
        }
    }
    if ( sampleBSDF ) {
        if ( ! light.isDeltaLight() ) {
            // return event.its->Ns;
            Spectrum f = bsdf->sampleF(event, uShading, false);
            scatteringPdf = event.pdf;
            if ( ! isBlack(f) && scatteringPdf != 0 ) {
                auto l = length(event.wi);
                vec3 worldShadowRayDir = event.toWorld(event.wi);
                Ray shaowRay(event.its->p, worldShadowRayDir);
                Spectrum Li = evalLightDirect(scene, light, shaowRay, medium, & lightPdf);
                if ( lightPdf == 0 ) return Ld;

                Float weight = 1;
                if ( ! isSpecualr(event.sampleType) ) {
                    weight = PowerHeuristic(scatteringPdf, lightPdf);
                    if ( ! sampleLgiht ) weight = 1;
                }

                if ( ! isBlack(Li) ) Ld += f * weight * Li / scatteringPdf;
                if ( hasNan(Ld) ) {
                    bsdf->f(event, false);
                    int k = 1;
                }
            }
        }
    }
    if ( hasNan(Ld) ) {

    }
    // assert(isBlack(Ld));
    return Ld;
}

Spectrum
Integrator::uniformSampleOneLight(SurfaceEvent & event,
                                  const Scene & scene,
                                  Sampler & sampler,
                                  const Distribution1D * lightDistrib,
                                  const Medium * medium) const {

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
    return estimateDirect(event, uScattering, * light, uLight,
                          scene, sampler, medium, false) / lightPdf;
}

Spectrum Integrator::uniformSampleAllLights(SurfaceEvent & event, const Scene & scene,
                                            Sampler & sampler, const Medium * medium) const {
    Spectrum L(0);
    for ( auto & light: scene.lights ) {
        L += estimateDirect(event, sampler.getNext2D(), * light, sampler.getNext2D(),
                            scene, sampler, medium, false);
    }
    return L;
}


SurfaceEvent makeLocalScatterEvent(const Intersection * its) {
    SurfaceEvent event;
    Frame frame = its->primitive->setTangentFrame(its);
    if( hasNan(frame.n)){
        its->primitive->setTangentFrame(its);
    }

    bool enableTwoSideShading = true;

    bool hitBackSide = dot(its->w, its->Ng) > 0;
    bool isTransmissive = its->bsdf->HasFlag(BSDF_TRANSMISSION);
    bool flippedFrame = false;
    //todo add this to config class
    if ( enableTwoSideShading && hitBackSide && ! isTransmissive ) {
        frame.n = - frame.n;
        frame.tangent = - frame.tangent;
        flippedFrame = true;
    }
    event.frame = frame;
    event.its = its;
    event.wo = event.toLocal(- its->w);
    event.flippedFrame = flippedFrame;
    return event;
}

std::unique_ptr < Distribution1D > Integrator::computeLightPowerDistrib(const Scene & scene) const {
    int nLights = scene.lights.size();
    Float * power = new Float[nLights];
    for ( int i = 0 ; i < nLights ; i ++ ) {
        power[i] = luminace(scene.lights[i]->Power());
    }
    return std::make_unique < Distribution1D >(power, nLights);
}

Spectrum
Integrator::evalLightDirect(const Scene & scene, const Light & light, Ray & ray,
                            const Medium * medium, Float * lightPdf) const {
    if ( light.isDeltaLight() )
        return Spectrum(0);
    std::optional < Intersection > lightIts = light.intersect(ray);
    if ( ! lightIts ) return Spectrum();
    Spectrum L;
    if ( light.flags & int(LightFlags::Area) ) L = lightIts->Le(- ray.d);
    if ( light.flags & int(LightFlags::Infinite) ) L = light.Le(ray);
    //avoid self shadow
    ray.farT -= lightIts->epsilon;
    if( isBlack(L)) return Spectrum(0);
    auto t = L * evalShadowDirect(scene, ray, medium);
    auto it = scene.intersect(ray);
    if ( hasNan(t) )
        int k = 1;
    if ( lightPdf ) * lightPdf = light.PdfLi(lightIts.value(), ray.o);
    return t;
}

Spectrum Integrator::evalShadowDirect(const Scene & scene, Ray ray, const Medium * medium) const {
    if ( ! medium ) {
        return scene.intersectP(ray) ? Spectrum(0) : Spectrum(1);
    }
    Spectrum Tr(1);
    std::optional < Intersection > its;
    Float tHit = ray.farT;
    while ( tHit > 0 ) {
        its = scene.intersect(ray);
        bool didHit = its.has_value();
        if ( didHit && ! its->bsdf->Pure(BSDF_FORWARD) )
            return Spectrum();
        Tr *= medium->TR(ray);
        medium = its->primitive->selectMedium(medium, dot(ray.d, its->Ng) > 0);
        tHit -= ray.farT;
        ray = Ray(its->p, ray.d, Constant::EPSILON, tHit);
    }
    return Tr;
}

Spectrum Integrator::volumeUniformSampleOneLight(VolumeEvent & event, const Medium * medium, const Scene & scene,
                                                 Sampler & sampler, const Distribution1D * lightDistrib) const {
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
    return volumeEstimateDirect(event, medium, uScattering, * light, uLight,
                                scene, sampler) / lightPdf;
}

Spectrum Integrator::volumeUniformSampleAllLights(VolumeEvent & event, const Medium * medium, const Scene & scene,
                                                  Sampler & sampler) const {
    return Spectrum();
}

Spectrum
Integrator::volumeEstimateDirect(VolumeEvent & event, const Medium * medium, const vec2 & uShading, const Light & light,
                                 const vec2 & uLight, const Scene & scene, Sampler & sampler) const {
    Spectrum Ld(0), Li;
    vec3 wi;
    Float lightPdf, distance, scatteringPdf;
    ///light sample
    // return event.p;
    Li = light.sampleLi(event.p, uLight, & wi, & lightPdf, & distance);
    Ray ray(event.p, wi, Constant::EPSILON, distance - Constant::EPSILON);
    // return Spectrum(distance);
    Spectrum p = event.phase->p(event.rayDir, wi);
    if ( ! isBlack(p) && ! isBlack(( Li = Li * evalShadowDirect(scene, ray, medium) )) ) {
        Float weight = 1;
        if ( ! light.isDeltaLight() ) weight = PowerHeuristic(lightPdf, scatteringPdf);
        Ld += Li * p * weight;
    }
    ///phase sample

    if ( ! light.isDeltaLight() ) {
        PhaseSample phaseSample;
        Spectrum f = event.phase->sampleP(event.rayDir, uShading, phaseSample);
        scatteringPdf = phaseSample.pdf;
        if ( ! isBlack(f) && scatteringPdf ) {
            Ray shaowRay(event.p, phaseSample.w, Constant::EPSILON);
            Li = evalLightDirect(scene, light, shaowRay, medium, & lightPdf);
            if ( lightPdf == 0 ) return Ld;
            Float weight = PowerHeuristic(scatteringPdf, lightPdf);
            if ( ! isBlack(Li) ) Ld += f * weight * Li / scatteringPdf;
        }
    }
    return Ld;
}

void SamplerIntegrator::render(const Scene & scene) const {
    //auto bmt =  std::make_shared <BitMapTexture<Spectrum>>("textures/envmap.hdr");
    //  bmt->LoadResources();
    auto uvbmt = std::make_shared < BitMapTexture < Spectrum>>("image18.png");
    uvbmt->LoadResources();

    const int tileSize = 16;
    ivec2 renderBounds = _camera->image->resoulation();
    int width = _camera->image->width();
    int height = _camera->image->height();
    ivec2 numTiles{( renderBounds.x + tileSize - 1 ) / tileSize, ( renderBounds.y + tileSize - 1 ) / tileSize};

    int num_threads = std::thread::hardware_concurrency();
    parallel_init(num_threads);

    int spp = scene.options.spp;
    int sppStep = scene.options.sppStep;

    ProgressReporter reporter(numTiles.x * numTiles.y);
    parallel_for([&](const vec2 & tile) {

        int x0 = tile[0] * tileSize;
        int x1 = std::min(x0 + tileSize, width);
        int y0 = tile[1] * tileSize;
        int y1 = std::min(y0 + tileSize, height);

        std::unique_ptr < Sampler > tileSampler = _sampler->clone();
        tileSampler->setSeed(tile.y * renderBounds.x + tile.x);
        for ( int y = y0 ; y < y1 ; y ++ ) {
            for ( int x = x0 ; x < x1 ; x ++ ) {
                for ( int s = 0 ; s < spp ; s ++ ) {
                    {
                        ///333 214  438 57 713 148
                        Ray ray = _camera->sampleRay(x, y, tileSampler->getNext2D());
                        Spectrum radiance = integrate(ray, scene, * tileSampler);
                        if ( hasNan(radiance) ) spdlog::error("Nan. Pixel {0} {1}", x, y);
                        _camera->image->addPixel(x, y, radiance);
                    }
                }
                _camera->image->dividePixel(x, y, spp);
            }
        }
        reporter.update(1);
    }, numTiles);
    _camera->image->postProgress();
    _camera->image->savePNG();

    parallel_cleanup();
}

bool russian(int depth, Float pdf,Sampler & sampler, vec3 &throguhPut) {

}

bool russian(int depth, Sampler &sampler, vec3 &beta) {
    float pdf = max(beta);
    if(depth<2 || pdf > 0.2)
        return false;
    if(sampler.getNext1D() < pdf )
    {
        beta/= pdf;
        return false;
    }
    return true;
}
