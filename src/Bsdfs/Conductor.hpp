#pragma once
#include "Reflection.hpp"
#include "MicrofacetDistribution.hpp"
class Conductor : public BSDF {
public:
    Conductor(vec3 eta,vec3 k) :
            BSDF(BXDFType(BSDF_SPECULAR | BSDF_REFLECTION)),
            m_eta(eta), m_k(k){}

    Conductor(std::string conductorName);
    Spectrum f(const SurfaceScatterEvent & event) const override;

    Float Pdf(const SurfaceScatterEvent & event) const override;

    Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;

    void LogInfo( ) const override;

protected:
    vec3 m_eta;
    vec3 m_k;
};


class RoughConductor : public BSDF{
public:
    RoughConductor(vec3 eta, vec3 k,  std::shared_ptr<MicrofacetDistribution> distrib,
                   std::shared_ptr < Texture<Float> > roughness,
                   std::shared_ptr < Texture<Float> > uroughness = nullptr,
                   std::shared_ptr < Texture<Float> > vroughness = nullptr
                   );
    void setCoundctorByName(const std::string & conductorName);
    Spectrum f(const SurfaceScatterEvent & event) const override;
    Float Pdf(const SurfaceScatterEvent & event) const override;
    Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;
    void LogInfo( ) const override;
    Float eta() const override;
protected:
    vec3 m_eta;
    vec3 m_k;
    std::shared_ptr<Texture<Float>> m_roughness,m_vRoughness, m_uRoughness;
    std::shared_ptr<MicrofacetDistribution> m_distrib;
};

