#include "Integrator.hpp"
#include "Bsdfs/Reflection.hpp"
#include "Common/ProgressReporter.h"
#include "Sampler/Sampler.hpp"
#include "Common/Parallel.h"
#include "Camera/Camera.hpp"
#include "PathIntegrator.hpp"
#include <thread>
#include <iostream>

static bool sampleBSDF = true;
static bool sampleLgiht = true;


//Integrator::Integrator(Json j) {
//
//}


Spectrum
Integrator::estimateDirect(SurfaceScatterEvent & event, const vec2 & uShading,
                           const Light & light, const vec2 & uLight,
                           const Scene & scene, Sampler & sampler, bool specular) const {

    BXDFType bsdfFlags =
            specular ? BSDF_ALL : BXDFType(BSDF_ALL & ~ BSDF_SPECULAR);
    event.sampleType = bsdfFlags;

    if(!event.its->bsdf->MatchesFlags(bsdfFlags))
        return Spectrum(0);

    vec3 wi;
    Float lightPdf = 0, scatteringPdf = 0;
    VisibilityTester visibility;

    Spectrum Ld(0);
    if ( sampleLgiht ) {
        Spectrum Li = light.sampleLi(* ( event.its ), uLight, & wi, & lightPdf, & visibility);

        if ( ! isBlack(Li) && lightPdf != 0  ) {
            if ( visibility.Unoccluded(scene)  )  {

                event.wi = event.toLocal(wi);
                scatteringPdf = event.its->bsdf->Pdf(event);
                Spectrum f = event.its->bsdf->f(event) * abs(event.wi.z);
                if ( ! isBlack(f) ) {
                    if ( light.isDeltaLight() ) Ld += f * Li;/// lightPdf;
                    else {
                        Float weight =
                                PowerHeuristic(lightPdf, scatteringPdf);
                        if ( ! sampleBSDF ) weight = 1;
                        Ld += f * weight * Li / lightPdf;
                        if( hasNan(Ld)){
                            int k=1;
                        }
                    }
                }
            }
            else {
            }
        }
    }

    if ( sampleBSDF )
        //sample bsdf
    {
        if ( ! light.isDeltaLight() ) {
           // return event.its->Ns;
            Spectrum f = event.its->bsdf->sampleF(event, uShading);
            f *= abs(event.wi.z);
            scatteringPdf = event.pdf;
            if ( !isBlack(f) && scatteringPdf )
            {
                Float weight = 1;
                vec3 worldShadowRayDir = event.toWorld(event.wi);
                Ray shaowRay(event.its->p, worldShadowRayDir,Constant::EPSILON);
                std::optional < Intersection > its = scene.intersect(shaowRay);
                bool isLightInfinite = ( light.flags & (int) LightFlags::Infinite );
                Spectrum  Li(0);
                if ( ( its.has_value() && its->primitive->areaLight.get() == & light ) || isLightInfinite ) {
                    if ( isLightInfinite ) {
                        if ( its.has_value() ) { return Ld; }
                        Intersection infinteIts;
                        its = std::make_optional(infinteIts);
                        its->w=-worldShadowRayDir;
                        Li = light.Le(shaowRay);
                    }
                    else {
                        Li = its->Le(-shaowRay.d);
                    }
                    if (!isSpecualr(event.sampleType) ) {
                        lightPdf = light.PdfLi(its.value(), event.its->p);
                        weight = PowerHeuristic(scatteringPdf, lightPdf);
                        if( isnan(weight))
                            lightPdf = light.PdfLi(its.value(), event.its->p);
                        if ( ! sampleLgiht ) weight = 1;
                    }
                } else {
                    return Ld;
                }

//                Spectrum Li = light.directLighting(its.value());
                Ld += f * weight * Li / scatteringPdf;
                if( hasNan(Ld)){
                    int k=1;
                }
            }
        }
    }
    if( hasNan(Ld)){

    }
   // assert(isBlack(Ld));
    return Ld;
}

Spectrum
Integrator::uniformSampleOneLight(SurfaceScatterEvent & event,
                                  const Scene & scene,
                                  Sampler & sampler,
                                  const Distribution1D * lightDistrib,
                                  bool handleMedia) const {

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
                          scene, sampler, handleMedia) / lightPdf;
}

Spectrum Integrator::uniformSampleAllLights(SurfaceScatterEvent & event, const Scene &scene,
                                            Sampler &sampler, bool handleMedia) const {
    Spectrum L(0);
    for(auto & light : scene.lights){
        L +=estimateDirect(event,sampler.getNext2D(), * light, sampler.getNext2D(),
                           scene, sampler, handleMedia);
    }
    return L;
}

SurfaceScatterEvent Integrator::makeLocalScatterEvent(const Intersection * its) const {
    SurfaceScatterEvent event;
    Frame frame = its->primitive->setTangentFrame(its);

    bool enableTwoSideShading = true ;

    bool hitBackSide = dot(its->w,its->Ng)>0;
    bool isTransmissive  = its->bsdf->hasFlag(BSDF_TRANSMISSION);
    bool flippedFrame = false;
    //todo add this to config class
    if(enableTwoSideShading && hitBackSide && !isTransmissive ){
        frame.n = -frame.n;
        frame.s = -frame.s;
        flippedFrame = true;
    }
    event.frame = frame;
    event.its = new Intersection(*its);
    event.wo =event.toLocal(-its->w);
    event.flippedFrame= flippedFrame;
    return event;
}

std::unique_ptr < Distribution1D > Integrator::computeLightPowerDistrib(const Scene & scene) const {
    int nLights = scene.lights.size();
    Float * power = new Float[nLights];
    for(int i=0;i<nLights;i++){
        power[i] = luminace(scene.lights[i]->Power());
    }
    return std::make_unique<Distribution1D>(power,nLights);
}

void SamplerIntegrator::render(const Scene & scene) const {
    const int tileSize = 16;
    ivec2 renderBounds = _camera->image->resoulation();
    int width = _camera->image->width();
    int height = _camera->image->height();
    ivec2 numTiles{(renderBounds.x+tileSize-1)/tileSize,(renderBounds.y+tileSize-1)/tileSize};

    int num_threads = std::thread::hardware_concurrency();
    parallel_init(num_threads);

    ProgressReporter reporter(numTiles.x * numTiles.y);
    parallel_for([&](const vec2 &tile){
        
        int x0 = tile[0] * tileSize;
        int x1 = std::min(x0 + tileSize,width);
        int y0 = tile[1] * tileSize;
        int y1 = std::min(y0 + tileSize, height);

        std::unique_ptr<Sampler> tileSampler = _sampler->clone();
        tileSampler->setSeed(tile.y * renderBounds.x + tile.x);
        int spp = scene.options.spp;
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                for (int s = 0; s < spp; s++) {
                    {
                        Ray ray = _camera->sampleRay(x,y,tileSampler->getNext2D());
                       // Spectrum  radiance = Spectrum(tileSampler->getNext1D(),tileSampler->getNext1D(),tileSampler->getNext1D());
                        Spectrum radiance = integrate(ray,scene,*tileSampler);
                        assert(!hasNan(radiance));
                        _camera->image->addPixel(x,y,radiance);
                    }
            }
                _camera->image->dividePixel(x,y,spp);
            }
    }
   // reporter.update();
    } ,numTiles);
    _camera->image->postProgress();
    _camera->image->savePNG();

    parallel_cleanup();
}
