#include "Reflection.hpp"

class Conductor : public BSDF {
public:
    Conductor(vec3 eta,vec3 k) :
            BSDF(BXDFType(BSDF_SPECULAR | BSDF_REFLECTION)),
            m_eta(eta), m_k(k){}

    Conductor(std::string conductor_name);
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
    Spectrum f(const SurfaceScatterEvent & event) const override;

    Float Pdf(const SurfaceScatterEvent & event) const override;

    Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;

    void LogInfo( ) const override;

    Float eta( ) const override;
};

