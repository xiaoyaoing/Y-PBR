#include "Reflection.hpp"
#include "Fresnel.hpp"
#include "ComplexIor.hpp"

#include "Common/Frame.hpp"
#include "Sampler/Warp.hpp"
#include "Texture/TextureFactory.hpp"

std::tuple<std::shared_ptr<Texture<Float>>, std::shared_ptr<Texture<Float>>, std::shared_ptr<Texture<Float>>>
loadRoughness(const Json& json) {
    std::shared_ptr<Texture<Float>> roughness = nullptr, uroughness = nullptr, vroughness = nullptr;
    roughness  = TextureFactory::LoadTexture<Float>(json, "roughness", 0.01f);
    uroughness = TextureFactory::LoadTexture<Float>(json, "urounghness");
    vroughness = TextureFactory::LoadTexture<Float>(json, "vrounghness");
    return std::make_tuple(roughness, uroughness, vroughness);
}

Spectrum LambertainR::f(const SurfaceEvent& event) const {
    if (event.wo.z < 0 || event.wi.z < 0) {
        return Spectrum();
    }
    Spectrum albedo = m_albedo->eval(event.its);
    return albedo * Constant::INV_PI * AbsCosTheta(event.wi);
}

Spectrum LambertainR::sampleF(SurfaceEvent& event, const vec2& u) const {

    if (event.wo.z < 0)
        return Spectrum(0);
    event.wi         = Warp::squareToCosineHemisphere(u);
    event.pdf        = Warp::squareToCosineHemispherePdf(event.wi);
    event.sampleType = BXDFType(BSDF_REFLECTION | BSDF_DIFFUSE);

    return m_albedo->eval(event.its) * Constant::INV_PI * AbsCosTheta(event.wi);
}

Spectrum LambertainT::f(const SurfaceEvent& event) const {
    return m_albedo->eval(event.its) * Constant::INV_PI;
}

Spectrum LambertainT::sampleF(SurfaceEvent& event, const vec2& u) const {
    return Spectrum();
}

Float LambertainR::Pdf(const SurfaceEvent& event) const {
    // if (event.wi.z <= 0.0f || event.wo.z <= 0.0f)
    //     return 0.0f;
    return Warp::squareToCosineHemispherePdf(event.wo);
}

Spectrum SpecularR::f(const SurfaceEvent& event) const {
    const vec3& wi = event.wi;
    const vec3& wo = event.wo;
    if (Frame::Reflect(wi) == wo)
        return Spectrum(1.0);
    return Spectrum(0.0);
}

Float SpecularR::Pdf(const SurfaceEvent& event) const {
    return 0;
}

Spectrum
SpecularR::sampleF(SurfaceEvent& event, const vec2& u) const {

    event.wi         = Frame::Reflect(event.wo);
    event.pdf        = 1;
    event.sampleType = BXDFType(BSDF_SPECULAR | BSDF_REFLECTION);
    return m_albedo->eval(event.its);
}

Spectrum BSDF::sampleF(SurfaceEvent& event, const vec2& u, bool adjoint) const {
    //        if (!MatchesFlags(event.requestType))
    //            return Spectrum();
    Spectrum fResult = sampleF(event, u);
    if (adjoint) {
        //        fResult *= std::abs(
        //                dot((event.toWorld(event.wo)), event.its->Ng) * event.wi.z /
        //                dot((event.toWorld(event.wi)), event.its->Ng) * event.wo.z);
        return fResult;

    } else {
        fResult *= sqr(eta(event));
    }

    if (hasNeg(fResult) || hasNan(fResult)) {
        DebugBreak();
    }
    return fResult;
}

Spectrum BSDF::f(const SurfaceEvent& event, bool adjoint) const {

    Spectrum fResult = f(event);
    if (adjoint) {
        return fResult;

    } else
        fResult *= sqr(eta(event));
    return fResult;
}