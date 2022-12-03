#include "PhaseFunction.hpp"
#include "Sampler/Warp.hpp"
Spectrum HenyeyGreenstein::p(const vec3 & wo, const vec3 & wi) const {
    return Spectrum();
}

Float HenyeyGreenstein::pdf(const vec3 & wo, const vec3 & wi) const {
    return 0;
}

Spectrum HenyeyGreenstein::sampleP(const vec3 & rayDir, const vec2 & u, PhaseSample & phaseSample) const {
    return Spectrum();
}

Spectrum IsotropicPhaseFunction::p(const vec3 & wo, const vec3 & wi) const {
    return Spectrum(Constant::INV_FOUR_PI);
}

Float IsotropicPhaseFunction::pdf(const vec3 & wo, const vec3 & wi) const {
    return Constant::INV_FOUR_PI;
}

Spectrum IsotropicPhaseFunction::sampleP(const vec3 & rayDir, const vec2 & u, PhaseSample & phaseSample) const {
    phaseSample.pdf = Constant::INV_FOUR_PI;
    phaseSample.w = Warp::squareToUniformSphere(u);
    return Spectrum(Constant::INV_FOUR_PI);
}
