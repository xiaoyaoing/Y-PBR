#pragma once

#include "Common/math.hpp"
#include "Colors/Spectrum.hpp"
#include "Common/Json.hpp"

struct PhaseSample {
    vec3  w;
    Float pdf;
};

class PhaseFunction {
public:
    virtual Spectrum p(const vec3& wo, const vec3& wi) const                                    = 0;
    virtual Float    pdf(const vec3& wo, const vec3& wi) const                                  = 0;
    virtual Spectrum sampleP(const vec3& rayDir, const vec2& u, PhaseSample& phaseSample) const = 0;
};

class IsotropicPhaseFunction : public PhaseFunction {
public:
    IsotropicPhaseFunction(const Json& json) {}
    Spectrum p(const vec3& wo, const vec3& wi) const override;

    Float pdf(const vec3& wo, const vec3& wi) const override;

    Spectrum sampleP(const vec3& rayDir, const vec2& u, PhaseSample& phaseSample) const override;
};

class HenyeyGreenstein : public PhaseFunction {
public:
    Spectrum p(const vec3& wo, const vec3& wi) const override;

    Float pdf(const vec3& wo, const vec3& wi) const override;

    Spectrum sampleP(const vec3& rayDir, const vec2& u, PhaseSample& phaseSample) const override;

protected:
    Float g;
};