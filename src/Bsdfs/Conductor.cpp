#include "Conductor.hpp"
#include "Fresnel.hpp"
#include "ComplexIor.hpp"

Conductor::Conductor(std::string conductor_name) : BSDF(BXDFType(BSDF_SPECULAR | BSDF_REFLECTION)){
    ComplexIorList::lookup(conductor_name,m_eta,m_k);
}

Float Conductor::Pdf(const SurfaceScatterEvent & event) const {
    return 0;
}

Spectrum Conductor::f(const SurfaceScatterEvent & event) const {
    return Spectrum();
}

Spectrum Conductor::sampleF(SurfaceScatterEvent & event, const vec2 & u) const {
    event.wi = Frame::Reflect(event.wo);
    event.pdf = 1;
    Spectrum  alebdo = m_albedo->Evaluate();
    Spectrum  f = alebdo  *  Fresnel::conductorReflectance(m_eta, m_k, event.wo.z);
    event.sampleType = m_type;
    return f;
}

void Conductor::LogInfo( ) const {
    spdlog::info("Specular Conductor albedo{0} eta{1} k{}", toColorStr(m_albedo->Evaluate()),
                 toColorStr(m_eta), toColorStr(m_k)
    );
}



Spectrum RoughConductor::f(const SurfaceScatterEvent & event) const {
    return Spectrum();
}

Float RoughConductor::Pdf(const SurfaceScatterEvent & event) const {
    return 0;
}

Spectrum
RoughConductor::sampleF(SurfaceScatterEvent & event, const vec2 & u) const {
    return Spectrum();
}

void RoughConductor::LogInfo( ) const {

}

Float RoughConductor::eta( ) const {
    return BSDF::eta();
}





