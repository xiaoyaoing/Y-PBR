#include "Reflection.hpp"
#include "MicrofacetDistribution.hpp"
class Plastic : public  BSDF{
public:
    Float Pdf(const SurfaceEvent & event) const override;

    Plastic(const std::shared_ptr < Texture < Spectrum>> & diffuseReflectance,
            const std::shared_ptr < Texture < Spectrum>> & specularReflectance, Float mIor):
            BSDF(BXDFType(BSDF_SPECULAR | BSDF_REFLECTION | BSDF_DIFFUSE)),
            diffuseReflectance(diffuseReflectance),specularReflectance(specularReflectance),m_ior(mIor)
    {}
protected:
    Spectrum sampleF(SurfaceEvent & event, const vec2 & u) const override;

    Spectrum f(const SurfaceEvent & event) const override;

    std::shared_ptr<Texture<Spectrum>> diffuseReflectance;
    std::shared_ptr<Texture<Spectrum>> specularReflectance;
    Float m_ior;
};

class Plastic1 : public BSDF{
public:
    Float Pdf(const SurfaceEvent & event) const override;

    Plastic1 (const Json & json);
protected:
    Spectrum sampleF(SurfaceEvent & event, const vec2 & u) const override;

    Spectrum f(const SurfaceEvent & event) const override;

     std::pair<Float,Float> specAndDiffuseProb(const SurfaceEvent & event) const;

    Float m_ior;
    Float _diffuseFresnel;
    Float _avgTransmittance;
    vec3 _scaledSigmaA;
};


class RoughPlastic : public  BSDF{
    std::shared_ptr<Texture<Spectrum>> diffuseReflectance;
    std::shared_ptr<Texture<Spectrum>> specularReflectance;
    std::shared_ptr<Texture<Float>> m_roughness,m_vRoughness, m_uRoughness;
    std::shared_ptr<MicrofacetDistribution> m_distrib;

    Float m_ior;
    Float _diffuseFresnel;
    Float _avgTransmittance;
    vec3 _scaledSigmaA;
public:
    RoughPlastic(const Json & json);
    RoughPlastic(const std::shared_ptr < Texture < Spectrum>> & diffuseReflectance,
                 const std::shared_ptr < Texture < Spectrum>> & specularReflectance, Float mIor,
                 const std::shared_ptr < MicrofacetDistribution > & mDistrib,
                 const std::shared_ptr < Texture < Float>> & mRoughness,
                 const std::shared_ptr < Texture < Float>> & mVRoughness,
                 const std::shared_ptr < Texture < Float>> & mURoughness);

    Spectrum f(const SurfaceEvent & event) const override;

    Float Pdf(const SurfaceEvent & event) const override;

    Spectrum sampleF(SurfaceEvent & event, const vec2 & u) const override;

    vec2 getAlphaXY(const SurfaceEvent & event) const;

};
