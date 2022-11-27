#include "Plastic.hpp"
#include "Fresnel.hpp"
#include "Sampler/Warp.hpp"
RoughPlastic::RoughPlastic(const std::shared_ptr < Texture < Spectrum>> & diffuseReflectance,
                           const std::shared_ptr < Texture < Spectrum>> & specularReflectance, Float mIor,
                           const std::shared_ptr < MicrofacetDistribution > & mDistrib,
                           const std::shared_ptr < Texture < Float>> & mRoughness,
                           const std::shared_ptr < Texture < Float>> & mVRoughness,
                           const std::shared_ptr < Texture < Float>> & mURoughness) :
                           BSDF(BXDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_DIFFUSE)),
                           diffuseReflectance(diffuseReflectance), specularReflectance(specularReflectance),
                           m_ior(mIor), m_distrib(mDistrib),
                           m_roughness(mRoughness),m_vRoughness(mVRoughness),m_uRoughness(mURoughness) {}
Spectrum RoughPlastic::f(const SurfaceScatterEvent & event) const {
    const vec3 & out = event.wo;
    const vec3 & in = event.wi;
    if(out.z <=0 || in.z<=0){
        return Spectrum(0);
    }
    vec3 wh = normalize((out + in));
    Float FOut = Fresnel::dielectricReflectance(1/m_ior,dot(out,wh));
    vec2 alphaXY = getAlphaXY(event);
    Float D = m_distrib->D(wh,alphaXY);
    Float G = m_distrib->G(event.wo,event.wi,alphaXY);
    Spectrum Ks = specularReflectance->Evaluate(event.its);
    Spectrum specularContrib = Ks * Spectrum(FOut * D * G / (4 * out.z));

    Float FIn = Fresnel::dielectricReflectance(1/m_ior,dot(in,wh));
    Spectrum Kd = diffuseReflectance->Evaluate(event.its);
    Spectrum diffuseContrib =Kd * (1-FOut) * (1-FIn) *(1/(m_ior * m_ior)) * * AbsCosTheta(event.wi) / Constant::PI;

    return specularContrib+diffuseContrib;
}

Float RoughPlastic::Pdf(const SurfaceScatterEvent & event) const {
    const vec3 & out = event.wo;
    const vec3 & in = event.wi;
    if( CosTheta(out) <=0 || CosTheta(in)<=0)
        return 0;
    vec3 wh = normalize(out + in);
    Spectrum Kd = diffuseReflectance->Evaluate(event.its);
    Spectrum Ks = specularReflectance->Evaluate(event.its);
    Float lS = luminace(Ks), lD = luminace(Kd);
    if (lS + lD <= 0) {
        return 0;
    }

    Float specProb = lS / (lS + lD);
    Float diffProb = 1 - specProb;
    vec2 alphaXY = getAlphaXY(event);
    Float D = m_distrib->D(wh,alphaXY);
    specProb *= D / (4 * absDot(out,wh));
    diffProb *= in.z / Constant::PI;
    return specProb+diffProb;
}

Spectrum RoughPlastic::sampleF(SurfaceScatterEvent & event, const vec2 & u) const {
    const vec3 & out = event.wo;
    if(CosTheta(out)<=0){
        event.pdf = 0;
        return {};
    }
    Spectrum Ks = specularReflectance->Evaluate(event.its);
    // We use the reflectance to choose between sampling the dielectric or diffuse layer.
    Spectrum Kd = diffuseReflectance->Evaluate(event.its);
    Float lS = luminace(Ks), lR = luminace(Kd);
    if (lS + lR <= 0) {
        event.pdf = 0;
        return {};
    }
    Float specProb = lS / (lS + lR);

    vec2 alphaxy = getAlphaXY(event);
    vec3 wh;
    if(u[0]<specProb){
        Float remapU0= (specProb-u[0])/specProb;
        vec2 newU(remapU0,u[1]);
        wh = m_distrib->Sample_wh(event.wo,newU,alphaxy);
        event.wi = Reflect(event.wo,wh);
        if(event.wi.z<=0)
            return Spectrum(0);
        event.sampleType= BXDFType(BSDF_REFLECTION | BSDF_GLOSSY);
        event.pdf = specProb *  m_distrib->D(wh,alphaxy) /(4 * absDot(out,wh)) + (1-specProb) * event.wi.z / Constant::PI;

    }
    else {
        Float remapU0= (u[0]-specProb)/(1-specProb);
        vec2 newU(remapU0,u[1]);
        event.wi = Warp::squareToCosineHemisphere(newU);
        if(event.wi.z<=0)
            return Spectrum(0);
        event.sampleType= BXDFType(BSDF_REFLECTION | BSDF_DIFFUSE);
        wh = normalize((event.wi+out));
        event.pdf = specProb *  m_distrib->D(wh,alphaxy) /(4 * absDot(out,wh)) + (1-specProb) * event.wi.z / Constant::PI;
    }

    Float FOut = Fresnel::dielectricReflectance(1/m_ior,dot(out,wh));
    Float D = m_distrib->D(wh,alphaxy);
    Float G = m_distrib->G(event.wo,event.wi,alphaxy);
    Spectrum specularContrib = Ks * Spectrum(FOut * D * G / (4 * out.z));
    Float FIn = Fresnel::dielectricReflectance(1/m_ior,dot(event.wi,wh));
    Spectrum diffuseContrib =Kd * (1-FOut) * (1-FIn) *(1/(m_ior * m_ior)) * AbsCosTheta(event.wi) / Constant::PI;
    return specularContrib+diffuseContrib;
}

void RoughPlastic::LogInfo( ) const {

}

vec2 RoughPlastic::getAlphaXY(const SurfaceScatterEvent & event) const {
    Float roughnessx = m_uRoughness? m_uRoughness->Evaluate(event.its) : m_roughness->Evaluate(event.its);
    Float roughnessy = m_vRoughness? m_vRoughness->Evaluate(event.its) : m_roughness->Evaluate(event.its);
    vec2 alphaXY = vec2(m_distrib->roughnessToAlpha(roughnessx),m_distrib->roughnessToAlpha(roughnessy));
    return alphaXY;
}
