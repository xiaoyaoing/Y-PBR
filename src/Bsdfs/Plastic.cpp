#include "Plastic.hpp"
#include "Fresnel.hpp"
#include "Sampler/Warp.hpp"
#include "Texture/TextureFactory.hpp"

RoughPlastic::RoughPlastic(const Json &json): BSDF(BXDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_DIFFUSE)) {
    m_ior = getOptional(json, "ior", 1.3);
    Float thickness = getOptional(json, "thickness", 1);
    Spectrum sigmaA = getOptional(json, "sigma_a", Spectrum(0.f));
    _scaledSigmaA = thickness * sigmaA;
    m_albedo = TextureFactory::LoadTexture<Spectrum>(json, "albedo", Spectrum(0.25));

    _avgTransmittance = std::exp(-2.0f * average(_scaledSigmaA));
    _diffuseFresnel = Fresnel::diffuseReflectance(m_ior, 10000);
    std::tie(m_roughness,m_uRoughness,m_vRoughness) = loadRoughness(json);
    m_distrib =LoadMicrofacetDistribution(getOptional(json, "distribution", std::string("beckmann")));
}



RoughPlastic::RoughPlastic(const std::shared_ptr<Texture<Spectrum>> &diffuseReflectance,
                           const std::shared_ptr<Texture<Spectrum>> &specularReflectance, Float mIor,
                           const std::shared_ptr<MicrofacetDistribution> &mDistrib,
                           const std::shared_ptr<Texture<Float>> &mRoughness,
                           const std::shared_ptr<Texture<Float>> &mVRoughness,
                           const std::shared_ptr<Texture<Float>> &mURoughness) :
        BSDF(BXDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_DIFFUSE)),
        diffuseReflectance(diffuseReflectance), specularReflectance(specularReflectance),
        m_ior(mIor), m_distrib(mDistrib),
        m_roughness(mRoughness), m_vRoughness(mVRoughness), m_uRoughness(mURoughness) {
//    Float thickness = getOptional, "thickness", 1);
//    Spectrum sigmaA = getOptional(json, "sigma_a", Spectrum(0.f));
    //_scaledSigmaA = thickness * sigmaA;


    _avgTransmittance = std::exp(-2.0f * average(_scaledSigmaA));
    _diffuseFresnel = Fresnel::diffuseReflectance(m_ior, 10000);
}

Spectrum RoughPlastic::f(const SurfaceEvent &event) const {
    const vec3 &out = event.wo;
    const vec3 &in = event.wi;
    if (out.z <= 0 || in.z <= 0) {
        return Spectrum(0);
    }
    vec3 wh = normalize((out + in));
    Float FOut = Fresnel::dielectricReflectance(1 / m_ior, dot(out, wh));
    vec2 alphaXY = getAlphaXY(event);
    Float D = m_distrib->D(wh, alphaXY);
    Float G = m_distrib->G(event.wo, event.wi, alphaXY);
    Spectrum specularContrib = Spectrum(FOut * D * G / (4 * out.z));

    Float FIn = Fresnel::dielectricReflectance(1 / m_ior, dot(in, wh));
    Spectrum albedo = m_albedo->eval(event.its);
    Spectrum diffuseContrib =
            albedo * (1 - FOut) * (1 - FIn) * (1 / (m_ior * m_ior)) * AbsCosTheta(event.wi) * Constant::INV_PI /(Spectrum(1.f)-_diffuseFresnel * albedo);
    if (max(_scaledSigmaA) > 0)
       diffuseContrib *= exp(_scaledSigmaA * (-1.0f / event.wo.z - 1.0f / event.wi.z));
    return specularContrib + diffuseContrib;
}

Float RoughPlastic::Pdf(const SurfaceEvent &event) const {
    
    const vec3 &out = event.wo;
    const vec3 &in = event.wi;
    if (CosTheta(out) <= 0 || CosTheta(in) <= 0)
        return 0;
    vec3 wh = normalize(out + in);
    auto F = Fresnel::dielectricReflectance(1/m_ior, event.wo.z);
    auto glossyProb = F;
    Float diffProb = (1 - glossyProb) * _avgTransmittance;
    glossyProb  /=(glossyProb + diffProb);
    diffProb = 1 - glossyProb;
    vec2 alphaXY = getAlphaXY(event);
    Float D = m_distrib->Pdf(event.wo, wh, alphaXY);
    glossyProb *= D / (4 * absDot(out, wh));
    diffProb *= in.z / Constant::PI;
    return glossyProb + diffProb;
}

Spectrum RoughPlastic::sampleF(SurfaceEvent &event, const vec2 &u) const {
    const vec3 &out = event.wo;
    if (CosTheta(out) <= 0) {
        event.pdf = 0;
        return {};
    }
  //  auto F = Fresnel::dielectricReflectance(1/m_ior, event.wo.z);
    Float FOut = Fresnel::dielectricReflectance(1 / m_ior, out.z);

    auto glossyProb = FOut;
    Float diffProb = (1 - glossyProb) * _avgTransmittance;
    glossyProb  /=(glossyProb + diffProb);
   

    vec2 alphaxy = getAlphaXY(event);
    vec3 wh;
    if (u[0] < glossyProb) {
        Float remapU0 = (glossyProb - u[0]) / glossyProb;
        vec2 newU(remapU0, u[1]);
        wh = m_distrib->Sample_wh(event.wo, newU, alphaxy);
        event.wi = Reflect(event.wo, wh);
        if (event.wi.z <= 0)
            return Spectrum(0);
        event.sampleType = BXDFType(BSDF_REFLECTION | BSDF_GLOSSY);
        event.pdf = glossyProb * m_distrib->D(wh, alphaxy) / (4 * absDot(out, wh)) +
                    (1 - glossyProb) * event.wi.z / Constant::PI;
    } else {
        Float remapU0 = (u[0] - glossyProb) / (1 - glossyProb);
        vec2 newU(remapU0, u[1]);
        event.wi = Warp::squareToCosineHemisphere(newU);
        if (event.wi.z <= 0)
            return Spectrum(0);
        event.sampleType = BXDFType(BSDF_REFLECTION | BSDF_DIFFUSE);
        wh = normalize((event.wi + out));
        event.pdf = glossyProb * m_distrib->D(wh, alphaxy) / (4 * absDot(out, wh)) +
                    (1 - glossyProb) * event.wi.z / Constant::PI;
    }
    auto albedo = m_albedo->eval(event.its);
    Float D = m_distrib->D(wh, alphaxy);
    Float G = m_distrib->G(event.wo, event.wi, alphaxy);
    Spectrum specularContrib =  Spectrum(FOut * D * G / (4 * out.z));
    Float FIn = Fresnel::dielectricReflectance(1 / m_ior, event.wi.z);
    Spectrum diffuseContrib =
            albedo * (1 - FOut) * (1 - FIn) * (1 / (m_ior * m_ior)) * AbsCosTheta(event.wi) /  (Spectrum(1)-_diffuseFresnel * albedo) / Constant::PI;
    if (max(_scaledSigmaA) > 0)
        diffuseContrib *= exp(_scaledSigmaA * (-1.0f / event.wo.z - 1.0f / event.wi.z));
    return specularContrib + diffuseContrib;
}

vec2 RoughPlastic::getAlphaXY(const SurfaceEvent &event) const {
    Float roughnessx = m_uRoughness ? m_uRoughness->eval(event.its) : m_roughness->eval(event.its);
    Float roughnessy = m_vRoughness ? m_vRoughness->eval(event.its) : m_roughness->eval(event.its);
    vec2 alphaXY = vec2(m_distrib->roughnessToAlpha(roughnessx), m_distrib->roughnessToAlpha(roughnessy));
    return alphaXY;
}

Float Plastic::Pdf(const SurfaceEvent &event) const {
    bool sampleSpecularR = (event.requestType & (BSDF_SPECULAR | BSDF_REFLECTION)) == (BSDF_SPECULAR | BSDF_REFLECTION);
    bool sampleDiffuseR = (event.requestType & (BSDF_DIFFUSE | BSDF_REFLECTION)) == (BSDF_DIFFUSE | BSDF_REFLECTION);
    Spectrum specularR = specularReflectance->eval(event.its), diffuseR = diffuseReflectance->eval(event.its);
    const vec3 &out = event.wo;
    const vec3 &in = event.wi;
    if (CosTheta(out) <= 0 || CosTheta(in) <= 0)
        return 0;
    Float lS = luminace(specularR), lD = luminace(diffuseR);
    if (lS + lD <= 0) {
        return 0;
    }
    Float specProb;
    if (sampleDiffuseR && sampleSpecularR) {
        specProb = lS / (lS + lD);
        specProb = Fresnel::dielectricReflectance(m_ior, event.wo.z);
        if (abs(out.z - in.z) < Constant::EPSILON) {
            return (1 - specProb) * Warp::squareToCosineHemispherePdf(in);
        } else
            return specProb;
    } else if (sampleSpecularR)
        return 1;
    else if (sampleDiffuseR)
        return Warp::squareToCosineHemispherePdf(in);
    return 0;
}

Spectrum Plastic::sampleF(SurfaceEvent &event, const vec2 &u) const {
    bool sampleSpecularR = (event.requestType & (BSDF_SPECULAR | BSDF_REFLECTION)) == (BSDF_SPECULAR | BSDF_REFLECTION);
    bool sampleDiffuseR = (event.requestType & (BSDF_DIFFUSE | BSDF_REFLECTION)) == (BSDF_DIFFUSE | BSDF_REFLECTION);

    Spectrum specularR = specularReflectance->eval(event.its), diffuseR = diffuseReflectance->eval(event.its);
    const vec3 &out = event.wo;
    if (CosTheta(out) <= 0) {
        event.pdf = 0;
        return {};
    }
    Float lS = luminace(specularR), lD = luminace(diffuseR);
    if (lS + lD <= 0) {
        event.pdf = 0;
        return Spectrum();
    }
    Float specProb;
    if (sampleDiffuseR && sampleSpecularR) {
        specProb = lS / (lS + lD);
        specProb = Fresnel::dielectricReflectance(m_ior, event.wo.z);
    } else if (sampleSpecularR)
        specProb = 1;
    else
        specProb = 0;
    if (u[0] < specProb) {
        event.wi = Frame::Reflect(out);
        event.sampleType = BXDFType(BSDF_REFLECTION | BSDF_SPECULAR);
        event.pdf = specProb;
        return specularR * Fresnel::dielectricReflectance(1 / m_ior, AbsCosTheta(out));
    } else {
        Float remapSample0 = (u[0] - specProb) / (1 - specProb);
        vec2 d(remapSample0, u[1]);
        event.wi = Warp::squareToCosineHemisphere(d);
        event.pdf = (1 - specProb) * Warp::squareToCosineHemispherePdf(event.wi);
        event.sampleType = BXDFType(BSDF_REFLECTION | BSDF_DIFFUSE);
        Float FOut = Fresnel::dielectricReflectance(1 / m_ior, AbsCosTheta(out));
        Float FIn = Fresnel::dielectricReflectance(1 / m_ior, AbsCosTheta(event.wi));
        return diffuseR * (1.f - FOut) * (1.f - FIn) * (1.f / (m_ior * m_ior)) * AbsCosTheta(event.wi) /
               Constant::PI;
    }
}

Spectrum Plastic::f(const SurfaceEvent &event) const {
    bool sampleSpecularR = (event.requestType & (BSDF_SPECULAR | BSDF_REFLECTION)) == (BSDF_SPECULAR | BSDF_REFLECTION);
    bool sampleDiffuseR = (event.requestType & (BSDF_DIFFUSE | BSDF_REFLECTION)) == (BSDF_DIFFUSE | BSDF_REFLECTION);
    Spectrum specularR = specularReflectance->eval(event.its), diffuseR = diffuseReflectance->eval(event.its);


    const vec3 &out = event.wo;
    const vec3 &in = event.wi;

    if (out.z <= 0 || in.z <= 0)
        return Spectrum(0);
    Float FOut = Fresnel::dielectricReflectance(1 / m_ior, out.z);

    if (abs(out.z - in.z) < Constant::EPSILON && sampleSpecularR)
        return specularR * FOut;
    else if (sampleDiffuseR) {

        Spectrum specularContrib(0);
        Float FIn = Fresnel::dielectricReflectance(1 / m_ior, in.z);
        return diffuseR * (1.f - FOut) * (1.f - FIn) * (1.f / (m_ior * m_ior)) * AbsCosTheta(in) / Constant::PI;
    }
}

Plastic1::Plastic1(const Json &json) : BSDF(BXDFType(BSDF_DIFFUSE | BSDF_REFLECTION | BSDF_SPECULAR)) {
    m_ior = getOptional(json, "ior", 1.3);
    Float thickness = getOptional(json, "thickness", 1);
    Spectrum sigmaA = getOptional(json, "sigma_a", Spectrum(0.f));
    _scaledSigmaA = thickness * sigmaA;
    m_albedo = TextureFactory::LoadTexture<Spectrum>(json, "albedo", Spectrum(0.25));

    _avgTransmittance = std::exp(-2.0f * average(_scaledSigmaA));
    _diffuseFresnel = Fresnel::diffuseReflectance(m_ior, 10000);
}

Float Plastic1::Pdf(const SurfaceEvent &event) const {
    bool sampleSpecularR = (event.requestType & (BSDF_SPECULAR | BSDF_REFLECTION)) == (BSDF_SPECULAR | BSDF_REFLECTION);
    bool sampleDiffuseR = (event.requestType & (BSDF_DIFFUSE | BSDF_REFLECTION)) == (BSDF_DIFFUSE | BSDF_REFLECTION);


    //  Spectrum specularR = specularReflectance->eval(event.its), diffuseR = diffuseReflectance->eval(event.its);
    const vec3 &out = event.wo;
    const vec3 &in = event.wi;
    if (CosTheta(out) <= 0 || CosTheta(in) <= 0)
        return 0;
    Float specProb, diffuseProb;
    if (sampleDiffuseR && sampleSpecularR) {
        specProb = Fresnel::dielectricReflectance(1/m_ior, out.z);
        diffuseProb = _avgTransmittance * (1 - specProb);
        specProb /= (specProb + diffuseProb);
        if (isMirrorReflect(out, in))
            return specProb;
        else return (1 - specProb) * Warp::squareToCosineHemispherePdf(in);
    } else if (sampleSpecularR) {
        return isMirrorReflect(out, in) ? 1 : 0;
    } else if (sampleDiffuseR)
        return Warp::squareToCosineHemispherePdf(in);
    else return 0;
}

Spectrum Plastic1::sampleF(SurfaceEvent &event, const vec2 &u) const {
    bool sampleSpecularR = (event.requestType & (BSDF_SPECULAR | BSDF_REFLECTION)) == (BSDF_SPECULAR | BSDF_REFLECTION);
    bool sampleDiffuseR = (event.requestType & (BSDF_DIFFUSE | BSDF_REFLECTION)) == (BSDF_DIFFUSE | BSDF_REFLECTION);

    const vec3 &out = event.wo;
    if (CosTheta(out) <= 0) {
        event.pdf = 0;
        return {};
    }

    Float specProb;
    auto F = Fresnel::dielectricReflectance(1/m_ior, out.z);
    if (sampleDiffuseR && sampleSpecularR) {
        specProb = F;
        Float diffuseProb = _avgTransmittance * (1 - F);
        specProb /= (specProb + diffuseProb);
    } else if (sampleSpecularR)
        specProb = 1;
    else
        specProb = 0;
    if (u[0] < specProb) {
        event.wi = Frame::Reflect(out);
        event.sampleType = BXDFType(BSDF_REFLECTION | BSDF_SPECULAR);
        event.pdf = specProb;
        return Spectrum(F);
    } else {
        Spectrum albedo = m_albedo->eval(event.its);
        Float remapSample0 = (u[0] - specProb) / (1 - specProb);
        vec2 d(remapSample0, u[1]);
        event.wi = Warp::squareToCosineHemisphere(d);
        event.pdf = (1-specProb) * Warp::squareToCosineHemispherePdf(event.wi);
        event.sampleType = BXDFType(BSDF_REFLECTION | BSDF_DIFFUSE);
        Float FOut = Fresnel::dielectricReflectance(1 / m_ior, AbsCosTheta(out));
        Float FIn = Fresnel::dielectricReflectance(1 / m_ior, AbsCosTheta(event.wi));
        auto brdf = albedo * (1.f - FOut) * (1.f - FIn) * (1.f / (m_ior * m_ior)) * event.wi.z * Constant::INV_PI /
                    (Spectrum(1.f) - albedo * _diffuseFresnel);
        if (max(_scaledSigmaA) > 0)
            brdf *= exp(_scaledSigmaA * (-1.0f / event.wo.z - 1.0f / event.wi.z));
        return brdf;
    }
}

Spectrum Plastic1::f(const SurfaceEvent &event) const {
    bool sampleSpecularR = (event.requestType & (BSDF_SPECULAR | BSDF_REFLECTION)) == (BSDF_SPECULAR | BSDF_REFLECTION);
    bool sampleDiffuseR = (event.requestType & (BSDF_DIFFUSE | BSDF_REFLECTION)) == (BSDF_DIFFUSE | BSDF_REFLECTION);
    auto albedo = m_albedo->eval(event.its);

    const vec3 &out = event.wo;
    const vec3 &in = event.wi;

    if (out.z <= 0 || in.z <= 0)
        return Spectrum(0);
    Float FOut = Fresnel::dielectricReflectance(1 / m_ior, out.z);

    if (isMirrorReflect(out, in) && sampleSpecularR)
        return Spectrum(FOut);
    else if (sampleDiffuseR) {
        Spectrum specularContrib(0);
        Float FIn = Fresnel::dielectricReflectance(1 / m_ior, in.z);
        auto brdf = albedo * (1.f - FOut) * (1.f - FIn) * (1.f / (m_ior * m_ior)) * event.wi.z * Constant::INV_PI /
                    (Spectrum(1.f) - albedo * _diffuseFresnel);
        if (max(_scaledSigmaA) > 0)
            brdf *= exp(_scaledSigmaA * (-1.0f / event.wo.z - 1.0f / event.wi.z));
        return brdf;
    }
}

std::pair<Float, Float> Plastic1::specAndDiffuseProb(const SurfaceEvent &event) const {
    auto specProb = Fresnel::dielectricReflectance(m_ior, event.wo.z);
    auto diffuseProb = (1 - specProb) * _avgTransmittance;
    specProb /= (specProb + diffuseProb);
    return {specProb, 1 - specProb};
}
