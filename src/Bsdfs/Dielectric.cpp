#include "Dielectric.hpp"


Spectrum Dielectric::f(const SurfaceScatterEvent & event) const {
    return Spectrum(0);
}

Float Dielectric::Pdf(const SurfaceScatterEvent & event) const {
    // avoid compute specular pdf
    return 0;
}

Spectrum Dielectric::sampleF(SurfaceScatterEvent & event, const vec2 & u) const {
    Spectrum  albedo = m_albedo->Evaluate();

    bool sampleT = enableT;
    Float  eta = event.wo.z < 0.0f ? ior : invIor;
    Float cosThetaT;
    Float F=Fresnel::DielectricReflectance(eta, std::abs(event.wo.z), cosThetaT);
    Float reflectionProbability=sampleT?F:1;
    if(u[0]<=reflectionProbability) {   //reflection
        event.wi =Frame::Reflect(event.wo);
        event.pdf = reflectionProbability;
        event.sampleType = BXDFType(BSDF_REFLECTION | BSDF_SPECULAR);
        return albedo * F / std::abs(event.wi.z);
    }
    else {
        if(reflectionProbability == 1.0){
            event.pdf =0;
            return  Spectrum(0);
        }
        event.wi = vec3(-eta * event.wo.x, -eta * event.wo.y, -std::copysign(cosThetaT, event.wo.z));

        event.pdf = 1-reflectionProbability;
        event.sampleType = BXDFType(BSDF_TRANSMISSION | BSDF_SPECULAR);
        return albedo * (1-F) / std::abs(event.wi.z);
    }
}

void Dielectric::LogInfo( ) const {
    spdlog::info("{DielectricBXDF ior:{0} albedo:{1} enableT{2}:",ior, toColorStr(m_albedo->Evaluate()),enableT);
}





Spectrum RoughDielectric::f(const SurfaceScatterEvent & event) const {
    return Spectrum();
}

Float RoughDielectric::Pdf(const SurfaceScatterEvent & event) const {
    return 0;
}

Spectrum
RoughDielectric::sampleF(SurfaceScatterEvent & event, const vec2 & u) const {
    return Spectrum();
}

void RoughDielectric::LogInfo( ) const {

}
