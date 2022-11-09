// MicrofacetDistribution Declarations
#pragma  once

#include "Common/math.hpp"
#include "Common/Json.hpp"








class MicrofacetDistribution {
public:
    // MicrofacetDistribution Public Methods
    virtual ~MicrofacetDistribution();
   // https://www.pbr-book.org/3ed-2018/Reflection_Models/Microfacet_Models
    virtual Float roughnessToAlpha(Float roughness) const =0;
    virtual Float D(const vec3 & wh, const vec2 & alphaxy) const = 0;
    virtual Float Lambda(const vec3 & w, const vec2 & alphaxy) const = 0;
    Float G1(const vec3 & w, const vec2 & alphaxy) const {
        //    if (Dot(w, wh) * CosTheta(w) < 0.) return 0.;
        return 1 / (1 + Lambda(w, alphaxy) );
    }
    virtual Float G(const vec3 & wo, const vec3 & wi, const vec2 & alphaxy) const {
        return 1 / ( 1 + Lambda(wo, alphaxy) + Lambda(wi, alphaxy) );
    }
    virtual vec3 Sample_wh(const vec3 & wo, const vec2 & u, const vec2 & alphaxy) const = 0;
    Float Pdf(const vec3 & wo, const vec3 & wh, const vec2 & alphaxy) const;
    virtual std::string ToString() const = 0;

protected:
    // MicrofacetDistribution Protected Methods
    MicrofacetDistribution(bool sampleVisibleArea)
            : sampleVisibleArea(sampleVisibleArea) {
    }

    // MicrofacetDistribution Protected Data
    const bool sampleVisibleArea;
};

class Beckmann : public  MicrofacetDistribution{
public:
    Beckmann(bool sampleVis = false): MicrofacetDistribution(sampleVis)
    {}
    Float roughnessToAlpha(float roughness) const override;
    Float D(const vec3 & wh, const vec2 & alphaxy) const override;
    vec3 Sample_wh(const vec3 & wo, const vec2 & u, const vec2 & alphaxy) const override;
    std::string ToString( ) const override;
protected:
    Float Lambda(const vec3 & w, const vec2 & alphaxy) const override;
};

class GGX : public MicrofacetDistribution {
public:
    GGX(): MicrofacetDistribution(true){}

    Float roughnessToAlpha(float roughness) const override;

    Float D(const vec3 & wh, const vec2 & alphaxy) const override;

    Float Lambda(const vec3 & w, const vec2 & alphaxy) const override;

    vec3 Sample_wh(const vec3 & wo, const vec2 & u, const vec2 & alphaxy) const override;

    std::string ToString( ) const override;
};

std::shared_ptr<MicrofacetDistribution>  LoadMicrofacetDistribution(const std::string & type);




//namespace  Mirofacet {
//    class Distribution{
//
//    };
//}
