#include "Light.hpp"
#include "Common/Json.hpp"

class DistantLight : public Light {
public:

    Float PdfLi(const Intersection & pShape, const vec3 & ref) const override;

    Spectrum
    sampleLi(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf, VisibilityTester * vis) const override;

    Spectrum Le(const Ray & ray) const override;

    Spectrum Power( ) override;

    void Preprocess(const Scene & scene) override;

    LightSampleResult sampleDirect(const vec2 & positionSample, const vec2 & dirSample) override;

    bool isDeltaLight( ) const override;

    DistantLight(const Json & json);

protected:
    Spectrum L;
    vec3 wLight;
    vec3 worldCenter;
    Float worldRadius;
};