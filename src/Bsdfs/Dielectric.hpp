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

    Spectrum f(const SurfaceEvent & event) const override;

    Float Pdf(const SurfaceEvent & event) const override;

    Spectrum sampleF(SurfaceEvent & event, const vec2 & u) const override;

    Float eta(const SurfaceEvent & event) const override {
        return 1;
        if(event.wi.z * event.wo.z >0)
            return 1;
        return event.wo.z<0?ior:invIor;
    }

protected:
    Float  ior,invIor;
    bool  enableT;
};


class RoughDielectric : public  BSDF{
public:
    RoughDielectric(Float _ior,  std::shared_ptr<MicrofacetDistribution> distrib ,
                    std::shared_ptr < Texture<Float> > roughness,
                    std::shared_ptr < Texture<Float> > uroughness = nullptr,
                    std::shared_ptr < Texture<Float> > vroughness = nullptr);

    Spectrum f(const SurfaceEvent & event) const override;
    Float Pdf(const SurfaceEvent & event) const override;

    Spectrum realF(const SurfaceEvent & event,bool reflect) const ;
    Float realPdf(const SurfaceEvent & event,bool reflect) const;

    Spectrum sampleF(SurfaceEvent & event, const vec2 & u) const override;

    vec2 getAlphaXY(const SurfaceEvent & event) const;

    Float eta(const SurfaceEvent & event) const override;

    Float m_ior;
    std::shared_ptr<MicrofacetDistribution> m_distrib;
    std::shared_ptr<Texture<Float>> m_roughness,m_vRoughness, m_uRoughness;

    //Dielectric dielectric;
};