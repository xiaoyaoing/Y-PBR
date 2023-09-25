#include "Dielectric.hpp"
#include "Fresnel.hpp"

Spectrum Dielectric::f(const SurfaceEvent & event) const {
    return Spectrum(0);
}

Float Dielectric::Pdf(const SurfaceEvent & event) const {
    // avoid compute specular pdf
    return 0;
}

Spectrum Dielectric::sampleF(SurfaceEvent & event, const vec2 & u) const {
    Spectrum  albedo = m_albedo->eval(event.its);

    bool sampleT = enableT;
    Float  eta = event.wo.z < 0.0f ? ior : invIor;
    Float cosThetaT;
    Float F=Fresnel::dielectricReflectance(eta, std::abs(event.wo.z), cosThetaT);
    Float reflectionProbability=sampleT?F:1;
    Spectrum f;
    if(u[0]<=reflectionProbability) {   //reflection
        event.wi =Frame::Reflect(event.wo);
        event.pdf = reflectionProbability;
        event.sampleType = BXDFType(BSDF_REFLECTION | BSDF_SPECULAR);
        f= albedo * F;
    }
    else {
        if(reflectionProbability == 1.0){
            event.pdf =0;
            return  Spectrum(0);
        }
        event.wi = vec3(-eta * event.wo.x, -eta * event.wo.y, -std::copysign(cosThetaT, event.wo.z));
        event.pdf = 1-reflectionProbability;
        event.sampleType = BXDFType(BSDF_TRANSMISSION | BSDF_SPECULAR);
        f = albedo * (1-F) ;
    }
    return f;
}

Float Dielectric::eta(const SurfaceEvent &event) const
{
        if(event.wi.z * event.wo.z >0)
            return 1;
        return event.wo.z<0?ior:invIor;
}


Spectrum RoughDielectric::f(const SurfaceEvent & event) const {
    const Spectrum albedo = m_albedo->eval(event.its->uv);
    const vec3 & out = event.wo;
    const vec3 & in = event.wi;
    bool reflect = out.z * in.z>0;
    Float eta = out.z >0 ? 1/m_ior : m_ior;
    vec3  wh;
    if (reflect) {
        wh = std::copysign(Float(1.0),out.z) * normalize(in + out);
    } else {
        wh = -normalize(in + out * eta);
    }
    Float F = Fresnel::dielectricReflectance(1/m_ior,dot(out,wh));

    vec2 alphaXY = getAlphaXY(event);

    Float D = m_distrib->D(wh,alphaXY);
    Float G = m_distrib->G(in,out,alphaXY);
    if (reflect) {
        return  albedo * F * D * G / ( 4* abs(out.z));
    }
    else {
        Float whDotIn =  dot(wh,in);
        Float whDotOut = dot(wh,out);
        Float sqrtDeom = eta * whDotOut  +  whDotIn;
        return  albedo * (1-F) *D * G  * std::abs(
                whDotIn * whDotOut  /
                (out.z * sqrtDeom * sqrtDeom)
        );
    }
}

Float RoughDielectric::Pdf(const SurfaceEvent & event) const {
    const vec3 & out = event.wo;
    const vec3 & in = event.wi;
    bool reflect = out.z * in.z > 0;
    Float eta = out.z >0 ? 1/m_ior : m_ior;
    vec3 wh;
    Float  pdf;
    if(reflect){
        wh = normalize(in+out) * Float(std::copysign(1,out.z));
    }
    else {
        wh = -normalize(in + out * eta);
    }
    Float  F =  Fresnel::dielectricReflectance(1/m_ior,dot(out,wh));
    vec2 alphaXY = getAlphaXY(event);
    Float whPdf = m_distrib->Pdf(out,wh,alphaXY);
    if(whPdf < 1e-50)
    {
        return 0;
    }
    if(reflect){
        pdf = F * whPdf / (4 * absDot(out,wh));
    }
    else {
        Float sqrtDenom = dot(out, wh) * eta +  dot(in, wh);
        Float dWhDWi =
                std::abs( dot(in, wh)) / (sqrtDenom * sqrtDenom);
        pdf =   whPdf * (1-F) * dWhDWi;
    }
    return pdf;
}

Spectrum RoughDielectric::sampleF(SurfaceEvent & event, const vec2 & u) const {
    const vec3  & out = event.wo;

    vec2 alphaXY  = getAlphaXY(event);

    vec3 wh = m_distrib->Sample_wh(out,u,alphaXY);
    Float whDotOut = dot(out, wh);
    Float  cosThetaT;
    Float  F = Fresnel::dielectricReflectance(1/m_ior,whDotOut,cosThetaT);
    Float  r = rand() % (10000 + 1) / (float)(10000 + 1);
    bool reflect = r<F;
    if(reflect){
        vec3 in = -out + 2 * dot(out, wh) * wh;
        event.wi = in;
        event.sampleType= BXDFType(BSDF_REFLECTION | BSDF_GLOSSY);
    }
    else {
        Float eta = whDotOut < 0.0f ? m_ior : 1.0f/m_ior;
        vec3 in = (eta * whDotOut - (whDotOut>0?1:-1)*cosThetaT)* wh - eta* out ;
        event.wi = in; 
        event.sampleType= BXDFType(BSDF_TRANSMISSION | BSDF_GLOSSY);
    }
    if(whDotOut<0){
        int k =1;
    }
    event.pdf= Pdf(event);
    return  f(event);
}

RoughDielectric::RoughDielectric(Float ior,  std::shared_ptr<MicrofacetDistribution> distrib,
                                 std::shared_ptr < Texture < Float>> roughness,
                                 std::shared_ptr < Texture < Float>> uroughness,
                                 std::shared_ptr < Texture < Float>> vroughness):
                                    m_distrib(distrib),
                                    BSDF(BXDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)),m_ior(ior),m_roughness(roughness),m_uRoughness(uroughness),m_vRoughness(vroughness)
                                 {
}

vec2 RoughDielectric::getAlphaXY(const SurfaceEvent & event) const {
    Float roughnessx = m_uRoughness ? m_uRoughness->eval(event.its) : m_roughness->eval(event.its);
    Float roughnessy = m_vRoughness ? m_vRoughness->eval(event.its) : m_roughness->eval(event.its);
    vec2 alphaXY = vec2(m_distrib->roughnessToAlpha(roughnessx),m_distrib->roughnessToAlpha(roughnessy));
    return alphaXY;
}

Float RoughDielectric::eta(const SurfaceEvent & event) const {
    if(event.wi.z * event.wo.z >0)
        return 1;
    return event.wo.z>0?1/m_ior:m_ior;
}
