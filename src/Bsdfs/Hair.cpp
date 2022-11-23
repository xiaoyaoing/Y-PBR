#include "Hair.h"
#include "Fresnel.hpp"
#include <iostream>

static Float trigInverse(Float x) {
    return sqrt(1 - x * x);
}

static Float I0(Float x) {
    float result = 1.0f;
    float xSq = x * x;
    float xi = xSq;
    float denom = 4.0f;
    for ( int i = 1 ; i <= 10 ; ++ i ) {
        result += xi / denom;
        xi *= xSq;
        denom *= 4.0f * float(( i + 1 ) * ( i + 1 ));
    }
    return result;
}


static float Phi(float gammaI, float gammaT, int p) {
    return 2.0f * p * gammaT - 2.0f * gammaI + p * Constant::PI;
}


static Float Gussian(Float beta, Float theta) {
    return std::exp(- theta * theta / ( 2.0f * beta * beta )) / ( std::sqrt(2.0f * Constant::PI) * beta );
}

static std::array < Spectrum, pMax + 1 > Ap(Float cosThetaO, Float eta, Float h,
                                            const Spectrum & T) {
    std::array < Spectrum, pMax + 1 > ap;
    // Compute $p=0$ attenuation at initial cylinder intersection
    Float cosGammaO = trigInverse(1 - h * h);
    Float cosTheta = cosThetaO * cosGammaO;
    Float f = Fresnel::dielectricReflectance(1 / eta, cosTheta);
    ap[0] = Spectrum(f);

    // Compute $p=1$ attenuation term
    ap[1] = sqr(1 - f) * T;

    // Compute attenuation terms up to $p=_pMax_$
    for ( int p = 2 ; p < pMax ; ++ p ) ap[p] = ap[p - 1] * T * f;

    // Compute attenuation term accounting for remaining orders of scattering
    ap[pMax] = ap[pMax - 1] * f * T / ( Spectrum(1.f) - T * f );
    return ap;
}

Float Hair::M(Float v, Float sinThetaO, Float sinThetaI, Float cosThetaO, Float cosThetaI) const {
    Float a = sinThetaO * sinThetaI / v;
    Float b = cosThetaI * cosThetaO / v;
    Float csch = 2 / ( exp(1 / v) - exp(- 1 / v) );
    return csch * exp(- b) * I0(- a);
}


Float Hair::NR(Float beta, Float cosO, Float phi, Float h) const {
    float gammaI = std::asin(clamp(h, - 1.0f, 1.0f));
    float deltaPhi = phi + 2.0f * gammaI;
    deltaPhi = std::fmod(deltaPhi, Constant::TWO_PI);
    if ( deltaPhi < 0.0f )
        deltaPhi += Constant::TWO_PI;
    return D(beta, deltaPhi) * Fresnel::dielectricReflectance(1.0f / _eta, cosO * cos(gammaI));
}

vec3 Hair::NP(Float beta, Float cosThetaD, Float phi, int p, Float h) const {
    float iorPrime = std::sqrt(_eta * _eta - ( 1.0f - cosThetaD * cosThetaD )) / cosThetaD;
    float cosThetaT = std::sqrt(1.0f - ( 1.0f - cosThetaD * cosThetaD ) * sqr(1.0f / _eta));
    vec3 sigmaAPrime = _sigmaA / cosThetaT;
    float gammaI = std::asin(clamp(h, - 1.0f, 1.0f));
    float gammaT = std::asin(clamp(h / iorPrime, - 1.0f, 1.0f));
    float l = 2.0f * std::cos(gammaT);

    float f = Fresnel::dielectricReflectance(1.0f / _eta, cosThetaD * trigInverse(h));
    vec3 T = exp(- sigmaAPrime * l);
    vec3 Aph = ( 1.0f - f ) * ( 1.0f - f ) * T;
    for ( int i = 1 ; i < p ; ++ i )
        Aph *= f * T;

    float deltaPhi = phi - Phi(gammaI, gammaT, p);
    deltaPhi = std::fmod(deltaPhi, Constant::TWO_PI);
    if ( deltaPhi < 0.0f )
        deltaPhi += Constant::TWO_PI;

    return Aph * D(beta, deltaPhi);
}

Spectrum Hair::f(const SurfaceScatterEvent & event) const {
    Float sinThetaO = event.wo.y;
    Float thetaO = std::asin(clamp(sinThetaO, - 1.0, 1.0));

    Float sinThetaI = event.wi.y;
    Float cosThetaI = sqrt(1 - sinThetaI * sinThetaI);
    Float thetaI = std::asin(clamp(sinThetaI, - 1.0, 1.0));

    Float thetaD = ( thetaI - thetaO ) * 0.5;
    Float cosThetaD = cos(thetaD);

    Float thetaOR = thetaO - 2.0 * _scaleAngle;
    Float thetaOTT = thetaO + _scaleAngle;
    Float thetaOTRT = thetaO + 4 * _scaleAngle;

    Float MR = M(_vR, sin(thetaOR), sinThetaI, cos(thetaOR), cosThetaI);
    Float MTT = M(_vTT, sin(thetaOTT), sinThetaI, cos(thetaOTT), cosThetaI);
    Float MTRT = M(_vTRT, sin(thetaOTRT), sinThetaI, cos(thetaOTRT), cosThetaI);

    Float phi = std::atan2(event.wi.x, event.wi.z);
    if ( phi < 0.0f )
        phi += Constant::TWO_PI;

    Float h = 0.5 + event.its->uv.y / 2;

    vec3 Nr = vec3(NR(_betaR, trigInverse(event.wo.y), phi, h));
    vec3 Ntt = NP(_betaTT, cosThetaD, phi, 1, h);
    vec3 Ntrt = NP(_betaTRT, cosThetaD, phi, 2, h);
    vec3 fsum = MR * Nr + MTT * Ntt + MTRT * Ntrt;
    return fsum / AbsCosTheta(event.wi);
}

Float Hair::Pdf(const SurfaceScatterEvent & event) const {
    return 0;
}

Spectrum Hair::sampleF(SurfaceScatterEvent & event, const vec2 & u) const {
    //Hair-samplineg requires 4 randoms.
    vec2 u0 = DemuxFloat(u[0]), u1 = DemuxFloat(u[1]);
    Float sinThetaO = event.wo.y;
    Float costhetaO = std::asin(clamp(sinThetaO, - 1.0, 1.0));
    Float thetaO = std::asin(clamp(sinThetaO, - 1, 1));
    //First choose a lobe to sample
    Float h = 0.5 + event.its->uv.x / 2;
    std::array < Float, pMax + 1 > apPdf = ComputeApPdf(costhetaO, h);
    int p;
    for ( p = 0 ; p < pMax ; ++ p ) {
        if ( u0[0] < apPdf[p] )
            break;
        u0[0] -= apPdf[p];
    }
    Float theta, v;
    if ( p == 0 ) {
        theta = thetaO - 2 * _scaleAngle;
        v = _betaR;
    } else if ( p == 1 ) {
        theta = thetaO + _scaleAngle;
        v = _betaTT;
    } else if ( p == 2 ) {
        theta = thetaO + 4 * _scaleAngle;
        v = _betaTRT;
    } else { throw ( "Invalid P" ); }

    Float sinThetaI = sampleM(v,sin(theta),cos(theta),u1[0],u1[1]);

    Float phiO = std::atan2(event.wo.x,event.wo.z);
    Float deltaphi;
    Float cosThetaD = cos(0.5 *(std::asin(sinThetaI)-thetaO));
    float iorPrime = std::sqrt(_eta * _eta - ( 1.0f - cosThetaD * cosThetaD )) / cosThetaD;
    float gammaI = std::asin(clamp(h, - 1.0f, 1.0f));
    float gammaT = std::asin(clamp(h / iorPrime, - 1.0f, 1.0f));

    return Spectrum();
}

void Hair::LogInfo( ) const {

}


Float Hair::D(Float beta, Float phi) const {
    float result = 0.0f;
    float delta;
    float shift = 0.0f;
    do {
        delta = Gussian(beta, phi + shift) + Gussian(beta, phi - shift - Constant::TWO_PI);
        result += delta;
        shift += Constant::TWO_PI;
    } while ( delta > 1e-4f );
    return result;
}

Hair::Hair(const Json & json) : BSDF(BXDFType(BSDF_GLOSSY | BSDF_TRANSMISSION | BSDF_REFLECTION)) {
    _scaleAngle = getOptional(json, "scale_angle", 2.5);
    _roughness = getOptional(json, "roughness", 0.3);
    Float melaninRatio = getOptional(json, "melanin_ratio", 1);
    Float melaninConcentration = getOptional(json, "melanin_concentration", 1.3);
    bool overrideSigmaA = containsAndGet(json, "sigma_a", _sigmaA);
    if ( ! overrideSigmaA ) {
        const vec3 eumelaninSigmaA = vec3(0.419f, 0.697f, 1.37f);
        const vec3 pheomelaninSigmaA = vec3(0.187f, 0.4f, 1.05f);

        _sigmaA = melaninConcentration * lerp(eumelaninSigmaA, pheomelaninSigmaA, melaninRatio);
    }

    _betaR = std::max(Constant::PI * 0.5f * _roughness, 0.04f);
    _betaTT = _betaR * 0.5f;
    _betaTRT = _betaR * 2.0f;

    _vR = _betaR * _betaR;
    _vTT = _betaTT * _betaTT;
    _vTRT = _betaTRT * _betaTRT;

    _scaleAngle = Angle::degToRad(_scaleAngle);


}

Float Hair::sampleM(Float v, Float sinThetaO, Float cosThetaO, Float xi1, Float xi2) const {
    float cosTheta = 1.0f + v * std::log(xi1 + ( 1.0f - xi1 ) * std::exp(- 2.0f / v));
    float sinTheta = trigInverse(cosTheta);
    float cosPhi = std::cos(Constant::TWO_PI * xi2);
    return - cosTheta * sinThetaO + sinTheta * cosPhi * cosThetaO;
}

std::array < Float, pMax + 1 > Hair::ComputeApPdf(Float cosThetaO, Float h) const {
    Float sinThetaO = trigInverse(cosThetaO);

    Float sinThetaT = sinThetaO / _eta;
    Float cosThetaT = trigInverse(sinThetaT);

    Float etap = std::sqrt(_eta * _eta - sqr(sinThetaO)) / cosThetaO;
    Float sinGammaT = h / etap;
    Float cosGammaT = trigInverse(cosGammaT);

    Spectrum T = exp(_sigmaA) * ( 2 * cosGammaT / cosThetaT );
    std::array < Spectrum, pMax + 1 > ap = Ap(cosThetaO, _eta, h, T);

    // Compute $A_p$ PDF from individual $A_p$ terms
    std::array < Float, pMax + 1 > apPdf;
    Float sumY =
            std::accumulate(ap.begin(), ap.end(), Float(0),
                            [](Float s, const Spectrum & ap) { return s + luminace(ap); });
    for ( int i = 0 ; i <= pMax ; ++ i ) apPdf[i] = luminace(ap[i]) / sumY;
    return apPdf;
}


