#include "Conductor.hpp"
#include "Fresnel.hpp"
#include "ComplexIor.hpp"

Conductor::Conductor(std::string conductorName) : BSDF(BXDFType(BSDF_SPECULAR | BSDF_REFLECTION)){
    ComplexIorList::lookup(conductorName, m_eta, m_k);
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
  //  f= Spectrum(1);
//    f= Spectrum(0.5);
    //f= Spectrum(0);
    event.sampleType = m_type;
    return f;
}

void Conductor::LogInfo( ) const {
    spdlog::info("Specular Conductor albedo{0} eta{1} k{}", toColorStr(m_albedo->Evaluate()),
                 toColorStr(m_eta), toColorStr(m_k)
    );
}



Spectrum RoughConductor::f(const SurfaceScatterEvent & event) const {
    if(!MatchesFlags(event.sampleType)){
        return Spectrum(0);
    }

    Float roughnessx = m_uRoughness? m_uRoughness->Evaluate(event.its) : m_roughness->Evaluate(event.its);
    Float roughnessy = m_vRoughness? m_vRoughness->Evaluate(event.its) : m_roughness->Evaluate(event.its);
    vec2 alphaxy = vec2(m_distrib->roughnessToAlpha(roughnessx),m_distrib->roughnessToAlpha(roughnessy));

    Float cosThetaO = AbsCosTheta(event.wo), cosThetaI = AbsCosTheta(event.wi);
    vec3 wh = event.wi + event.wo;
    // Handle degenerate cases for microfacet reflection
    if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0.);
    if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.);
    wh = normalize(wh);
    // For the Fresnel call, make sure that wh is in the same hemisphere
    // as the surface normal, so that TIR is handled correctly.
    Float cosI = dot(event.wi,faceForward(wh,vec3(0,0,1)));
    Spectrum F= Fresnel::conductorReflectance(m_eta,m_k,cosI);
    return  m_albedo->Evaluate(event.its) * m_distrib->D(wh,alphaxy)  * F * m_distrib->G(event.wo, event.wi,alphaxy)/ (4 * cosThetaI * cosThetaO);
}

Float RoughConductor::Pdf(const SurfaceScatterEvent & event) const {
    if (!SameHemisphere(event.wo,event.wi)) return 0;
    Float roughnessx = m_uRoughness? m_uRoughness->Evaluate(event.its) : m_roughness->Evaluate(event.its);
    Float roughnessy = m_vRoughness? m_vRoughness->Evaluate(event.its) : m_roughness->Evaluate(event.its);
    vec2 alphaxy = vec2(m_distrib->roughnessToAlpha(roughnessx),m_distrib->roughnessToAlpha(roughnessy));

    vec3 wh = normalize(event.wo + event.wi);
    return m_distrib->Pdf(event.wo, wh,alphaxy) / (4 * dot(event.wo, wh));
}

Spectrum
RoughConductor::sampleF(SurfaceScatterEvent & event, const vec2 & u) const {
    if(!MatchesFlags(event.sampleType)){
        return Spectrum(0);
    }
    Float roughnessx = m_uRoughness? m_uRoughness->Evaluate(event.its) : m_roughness->Evaluate(event.its);
    Float roughnessy = m_vRoughness? m_vRoughness->Evaluate(event.its) : m_roughness->Evaluate(event.its);

    vec2 alphaxy = vec2(m_distrib->roughnessToAlpha(roughnessx),m_distrib->roughnessToAlpha(roughnessy));
    vec3 wh = m_distrib->Sample_wh(event.wo,u,alphaxy);
    event.wi = Reflect(event.wo,wh);
    if (!SameHemisphere(event.wi, event.wo))
        return Spectrum(0.f);

    event.sampleType = m_type;
    // Compute PDF of _wi_ for microfacet reflection
    //这里除法是要转换分布（法线到入射光）
    event.pdf =m_distrib->Pdf(event.wo, wh,alphaxy) / (4 * dot(event.wo, wh));
    return  f(event);

}

void RoughConductor::LogInfo( ) const {

}

Float RoughConductor::eta( ) const {
    return BSDF::eta();
}

RoughConductor::RoughConductor(vec3 eta, vec3 k, MircroDistributionEnum distribType,
                               std::shared_ptr < Texture<Float> > roughness, std::shared_ptr < Texture<Float> > uroughness,
                               std::shared_ptr < Texture<Float> > vroughness):
                               BSDF(BXDFType(BSDF_GLOSSY | BSDF_REFLECTION)),
                               m_eta(eta), m_k(k),m_roughness(roughness),m_uRoughness(uroughness),m_vRoughness(vroughness)
                               {
    if(distribType ==  GGx){

    }
    else if(distribType == Beckmann){
        this->m_distrib = std::make_shared<BeckmannDistribution>();
    }

}

void RoughConductor::setCoundctorByName(const std::string & conductorName) {
    ComplexIorList::lookup(conductorName, m_eta, m_k);
}




