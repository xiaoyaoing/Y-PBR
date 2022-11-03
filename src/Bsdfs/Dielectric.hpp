#pragma  once

#include "Reflection.hpp"
#include "MicrofacetDistribution.hpp"

class Dielectric : public  BSDF {
public:

    Dielectric(Float ior,bool enableT=true) : BSDF(BXDFType(BSDF_SPECULAR | BSDF_REFLECTION |
                                                            (enableT?BSDF_TRANSMISSION:0)) ),
                                              ior(ior), enableT(enableT){
        invIor = 1.0f/ior;
    }

    Spectrum f(const SurfaceScatterEvent & event) const override;

    Float Pdf(const SurfaceScatterEvent & event) const override;

    Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;

    void LogInfo( ) const override;

    Float eta( ) const override {
        return ior;
    }

private:
    Float  ior,invIor;
    bool  enableT;
};


class RoughDielectric : public  BSDF{
public:
    Spectrum f(const SurfaceScatterEvent & event) const override;

    Float Pdf(const SurfaceScatterEvent & event) const override;

    Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;

    void LogInfo( ) const override;
};