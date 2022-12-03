#include "Reflection.hpp"
#include "Fresnel.hpp"
#include "ComplexIor.hpp"

#include "Common/Frame.hpp"
#include "Sampler/Warp.hpp"

#include <spdlog/spdlog.h>

Spectrum LambertainR::f(const SurfaceEvent & event) const {
    if ( event.wo.z < 0 || event.wi.z < 0 ) {
        return Spectrum();
    }
    Spectrum albedo = m_albedo->Evaluate(event.its);
    return albedo * Constant::INV_PI * AbsCosTheta(event.wi);
}


void LambertainR::LogInfo( ) const {
    Spectrum albedo = m_albedo->Evaluate();
    spdlog::info("{0} albedo {1} {2} {3}", "LambertainReflection", albedo.x, albedo.y, albedo.z);
}


Spectrum LambertainR::sampleF(SurfaceEvent & event, const vec2 & u) const {

    event.wi = Warp::squareToUniformHemisphere(u);
    event.pdf = Warp::squareToUniformHemispherePdf(event.wi);
    event.sampleType = BXDFType(BSDF_REFLECTION | BSDF_DIFFUSE);

    return m_albedo->Evaluate(event.its) * Constant::INV_PI * AbsCosTheta(event.wi);
}


Spectrum LambertainT::f(const SurfaceEvent & event) const {
    return m_albedo->Evaluate(event.its) * Constant::INV_PI;
}

void LambertainT::LogInfo( ) const {

}


Spectrum LambertainT::sampleF(SurfaceEvent & event, const vec2 & u) const {
    return Spectrum();
}


Float LambertainR::Pdf(const SurfaceEvent & event) const {
    if ( event.wi.z <= 0.0f || event.wo.z <= 0.0f )
        return 0.0f;
    return Warp::squareToUniformHemispherePdf(event.wo);
}

Spectrum SpecularR::f(const SurfaceEvent & event) const {
    const vec3 & wi = event.wi;
    const vec3 & wo = event.wo;
    if ( Frame::Reflect(wi) == wo )
        return Spectrum(1.0);
    return Spectrum(0.0);
}

Float SpecularR::Pdf(const SurfaceEvent & event) const {
    return 0;
}


void SpecularR::LogInfo( ) const {
    spdlog::info("Specular Reflection");
}


Spectrum
SpecularR::sampleF(SurfaceEvent & event, const vec2 & u) const {

    event.wi = Frame::Reflect(event.wo);
    event.pdf =  1;
    event.sampleType = BXDFType(BSDF_SPECULAR | BSDF_REFLECTION);
    return  m_albedo->Evaluate(event.its);
}

