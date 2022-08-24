// MicrofacetDistribution Declarations
#include "Common/math.hpp"

class MicrofacetDistribution {
public:
    // MicrofacetDistribution Public Methods
    virtual ~MicrofacetDistribution();
    virtual Float D(const vec3 &wh) const = 0;
    virtual Float Lambda(const vec3 &w) const = 0;
    Float G1(const vec3 &w) const {
        //    if (Dot(w, wh) * CosTheta(w) < 0.) return 0.;
        return 1 / (1 + Lambda(w));
    }
    virtual Float G(const vec3 &wo, const vec3 &wi) const {
        return 1 / (1 + Lambda(wo) + Lambda(wi));
    }
    virtual vec3 Sample_wh(const vec3 &wo, const vec2 &u) const = 0;
    Float Pdf(const vec3 &wo, const vec3 &wh) const;
    virtual std::string ToString() const = 0;

protected:
    // MicrofacetDistribution Protected Methods
    MicrofacetDistribution(bool sampleVisibleArea)
            : sampleVisibleArea(sampleVisibleArea) {}

    // MicrofacetDistribution Protected Data
    const bool sampleVisibleArea;
};