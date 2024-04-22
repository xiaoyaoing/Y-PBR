#include "VolPathIntegrator.h"
#include "Bsdfs/Reflection.hpp"
#include "Mediums/Medium.hpp"
#include "TraceHelper.h"

vec3 VolPathIntegrator::integrate(const Ray& ray, const Scene& scene, Sampler& sampler) const {
    Spectrum      L(0);
    const Medium* medium;
    medium = _camera->_medium.get();
    VolumeEvent                 voulumeEvent{};
    SurfaceEvent                surfaceEvent;
    bool                        specularBounce = true;
    int                         bounce;
    Spectrum                    beta(1);
    Ray                         _ray(ray);
    std::optional<Intersection> temp;
    bool                        gemoback;
    bool                        tempHitSurface;
    std::optional<Intersection> its;
    for (bounce = 0; bounce < maxBounces; bounce++) {
        its                    = scene.intersect(_ray);
        bool foundIntersection = its.has_value();
        bool hitSurface        = true;
        if (!medium && !foundIntersection)
            break;

        if (medium) {
            medium->sampleDistance(_ray, sampler, voulumeEvent);
            beta *= medium->sampleDistance(_ray, sampler, voulumeEvent);
            if (isBlack(beta)) break;
            hitSurface = voulumeEvent.exited;
            if (!foundIntersection && hitSurface)
                break;
        }

        if (hitSurface) {

            if (specularBounce) {
                if (its.has_value())
                    L += beta * its->Le(-_ray.d);
            }
            surfaceEvent = makeLocalScatterEvent(&(its.value()));
            if (its->bsdf->Pure(BSDF_FORWARD)) {
                _ray = surfaceEvent.sctterRay(_ray.d);
            } else {
                if (its->bsdf->MatchesFlags(BXDFType(BSDF_NO_SPECULAR)) && bounce < maxBounces - 1) {
                    auto l = uniformSampleAllLights(surfaceEvent, scene, sampler, medium);
                    L += beta * l;
                }
                surfaceEvent.requestType = BSDF_ALL;
                Spectrum f               = its->bsdf->sampleF(surfaceEvent, sampler.getNext2D(), false);
                if (isBlack(f) || surfaceEvent.pdf == 0)
                    break;
                BXDFType flags = surfaceEvent.sampleType;
                specularBounce = (flags & BSDF_SPECULAR) != 0;
                if (specularBounce) {
                    DebugBreak();
                }
                beta *= f / surfaceEvent.pdf;
                _ray = surfaceEvent.sctterRay(surfaceEvent.toWorld(surfaceEvent.wi));
            }
            medium = its->primitive->selectMedium(medium, dot(_ray.d, its->Ng) < 0);
        } else {
            specularBounce = false;
            auto l         = volumeUniformSampleOneLight(voulumeEvent, medium, scene, sampler, lightDistribution.get());
            L += beta *
                 l;
            // return beta * l;
            if ((luminace(L) > 2)) {
                DebugBreak();
            }

            PhaseSample phaseSample;
            Spectrum    p = voulumeEvent.phase->sampleP(_ray.d, sampler.getNext2D(), phaseSample);
            beta *= p / phaseSample.pdf;
            _ray = Ray(voulumeEvent.p, phaseSample.w);
        }
        //russian prob
        Float roulettePdf = std::max(beta.x, std::max(beta.y, beta.z));
        if (bounce > 2 && roulettePdf < 0.1) {
            if (sampler.getNext1D() < roulettePdf)
                beta /= roulettePdf;
            else
                break;
        }
        if (bounce == 2 && luminace(L) > 2) {
            DebugBreak();
        }
    }

    if (specularBounce)
        for (auto light : scene.lights) {
            if (light->flags == int(LightFlags::Infinite)) {
                L += beta * light->Le(_ray);
            }
        }
    if (luminace(L) == 0) {
        DebugBreak();
    }
    if (bounce == 3) {
    }
    return L;
}