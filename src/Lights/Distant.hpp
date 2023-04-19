#include "Light.hpp"
#include "Common/Json.hpp"
#include "PositionAndDirectionSample.h"

class DistantLight : public Infinite {
public:

    Float PdfLi(const Intersection & pShape, const vec3 & ref) const override;

    Spectrum
    sampleLi(const vec3 & ref, const vec2 & u, vec3 * wi, Float * pdf, Float * distance) const override;

    Spectrum Le(const Ray & ray) const override;

    Spectrum Power( ) override;


    PositionAndDirectionSample sampleDirect(const vec2 & positionSample, const vec2 & dirSample) override;

    bool isDeltaLight( ) const override;

    DistantLight(const Json & json);

protected:
    Spectrum L;
    vec3 wLight;

};