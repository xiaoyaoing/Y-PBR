#include "Reflection.hpp"
#include "MicrofacetDistribution.hpp"
class Plastic : public  BSDF{

};


class RoughPlastic : public  BSDF{
    std::shared_ptr<Texture<Spectrum>> diffuseReflectance;
    std::shared_ptr<Texture<Spectrum>> specularReflectance;
    Float m_ior;
    std::shared_ptr<Texture<Float>> m_roughness,m_vRoughness, m_uRoughness;
    std::shared_ptr<MicrofacetDistribution> m_distrib;
public:
    RoughPlastic(const std::shared_ptr < Texture < Spectrum>> & diffuseReflectance,
                 const std::shared_ptr < Texture < Spectrum>> & specularReflectance, Float mIor,
                 const std::shared_ptr < MicrofacetDistribution > & mDistrib,
                 const std::shared_ptr < Texture < Float>> & mRoughness,
                 const std::shared_ptr < Texture < Float>> & mVRoughness,
                 const std::shared_ptr < Texture < Float>> & mURoughness);

    Spectrum f(const SurfaceScatterEvent & event) const override;

    Float Pdf(const SurfaceScatterEvent & event) const override;

    Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;

    void LogInfo( ) const override;

    vec2 getAlphaXY(const SurfaceScatterEvent & event) const;

};
