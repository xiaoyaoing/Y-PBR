#include "Hair.h"
#include "Fresnel.hpp"
#include <iostream>


bool useGussianAmz = true;


static Float trigInverse(Float x) {
    return std::min(std::sqrt(std::max(1.0f - x * x, 0.0f)), 1.0f);

}

static Float I0(Float x) {
    Float result = 1.0f;
    Float xSq = x * x;
    Float xi = xSq;
    Float denom = 4.0f;
    for (int i = 1; i <= 10; ++i) {
        result += xi / denom;
        xi *= xSq;
        denom *= 4.0f * Float((i + 1) * (i + 1));
    }
    return result;
}

static Float getH(const SurfaceEvent &event) {
    //return std::fabs(event.its->uv.y);
    return 2 * std::fabs(event.its->uv.y) - 1;
}


static Float Phi(Float gammaI, Float gammaT, int p) {
    return 2.0f * p * gammaT - 2.0f * gammaI + p * Constant::PI;
}


static Float Gussian(Float beta, Float theta) {
    return std::exp(-theta * theta / (2.0f * beta * beta)) / (std::sqrt(2.0f * Constant::PI) * beta);
}


inline Float Logistic(Float x, Float s) {
    x = std::abs(x);
    return std::exp(-x / s) / (s * sqr(1 + std::exp(-x / s)));
}

inline Float LogisticCDF(Float x, Float s) {
    return 1 / (1 + std::exp(-x / s));
}


inline Float TrimmedLogistic(Float x, Float s, Float a, Float b) {
    while (x > Constant::PI) x -= 2 * Constant::PI;
    while (x < -Constant::PI) x += 2 * Constant::PI;
    return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
}

static Float SampleTrimmedLogistic(Float u, Float s, Float a, Float b) {
    Float k = LogisticCDF(b, s) - LogisticCDF(a, s);
    Float x = -s * std::log(1 / (u * k + LogisticCDF(a, s)) - 1);
    return clamp(x, a, b);
}

static std::array<Spectrum, pMax + 1> Ap(Float cosThetaO, Float eta, Float h,
                                         const Spectrum &T) {
    std::array<Spectrum, pMax + 1> ap;
    // Compute $p=0$ attenuation at initial cylinder intersection
    Float cosGammaO = trigInverse(1 - h * h);
    Float cosTheta = cosThetaO * cosGammaO;
    Float f = Fresnel::dielectricReflectance(1 / eta, cosTheta);
    ap[0] = Spectrum(f);
    // Compute $p=1$ attenuation term
    ap[1] = sqr(1 - f) * T;
    // Compute attenuation terms up to $p=_pMax_$
    for (int p = 2; p < pMax; ++p) ap[p] = ap[p - 1] * T * f;
    // Compute attenuation term accounting for remaining orders of scattering
    // ap[pMax] = ap[pMax - 1] * f * T / ( Spectrum(1.f) - T * f );
    return ap;
}

static float logI0(float x) {
    if (x > 12.0f)
        // More stable evaluation of log(I0(x))
        // See also https://publons.com/discussion/12/
        return x + 0.5f * (std::log(1.0f / (Constant::TWO_PI * x)) + 1.0f / (8.0f * x));
    else
        return std::log(I0(x));
}

Float Hair::M(Float v, Float sinThetaO, Float sinThetaI, Float cosThetaO, Float cosThetaI) const {
    Float b = sinThetaO * sinThetaI / v;
    Float a = cosThetaI * cosThetaO / v;
    if (v < 0.1)
        return std::exp(-b + logI0(a) - 1.0f / v + 0.6931f + std::log(1.0f / (2.0f * v)));
    Float csch = 2 / (exp(1 / v) - exp(-1 / v));
    return csch * exp(-b) * I0(a);
}

Float Hair::NR(Float beta, Float cosO, Float phi, Float h) const {
    Float gammaI = std::asin(clamp(h, -1.0f, 1.0f));
    Float deltaPhi = phi + 2.0f * gammaI;
    deltaPhi = std::fmod(deltaPhi, Constant::TWO_PI);
    if (deltaPhi < 0.0f)
        deltaPhi += Constant::TWO_PI;
    // return D(beta, deltaPhi) * Fresnel::dielectricReflectance(1.0f / _eta, cosO * cos(gammaI));
    return TrimmedLogistic(deltaPhi, beta, -Constant::PI, Constant::PI) *
           Fresnel::dielectricReflectance(1.0f / _eta, cosO * cos(gammaI));
}

vec3 Hair::NP(Float beta, Float cosThetaD, Float phi, int p, Float h) const {
    Float iorPrime = std::sqrt(_eta * _eta - (1.0f - cosThetaD * cosThetaD)) / cosThetaD;
    Float cosThetaT = std::sqrt(1.0f - (1.0f - cosThetaD * cosThetaD) * sqr(1.0f / _eta));
    vec3 sigmaAPrime = _sigmaA / cosThetaT;
    Float gammaI = std::asin(clamp(h, -1.0f, 1.0f));
    Float gammaT = std::asin(clamp(h / iorPrime, -1.0f, 1.0f));
    Float l = 2.0f * std::cos(gammaT);

    Float f = Fresnel::dielectricReflectance(1.0f / _eta, cosThetaD * trigInverse(h));
    vec3 T = exp(-sigmaAPrime * l);
    vec3 Aph = (1.0f - f) * (1.0f - f) * T;
    for (int i = 1; i < p; ++i)
        Aph *= f * T;

    Float deltaPhi = phi - Phi(gammaI, gammaT, p);
    deltaPhi = std::fmod(deltaPhi, Constant::TWO_PI);
    if (deltaPhi < 0.0f)
        deltaPhi += Constant::TWO_PI;
    // return Aph * D(beta, deltaPhi);
    return Aph * vec3(TrimmedLogistic(deltaPhi, beta, -Constant::PI, Constant::PI));
}

Spectrum Hair::f(const SurfaceEvent &event) const {
    Float sinThetaO = event.wo.y;
    Float thetaO = std::asin(clamp(sinThetaO, -1.0, 1.0));

    Float sinThetaI = event.wi.y;
    Float cosThetaI = sqrt(1 - sinThetaI * sinThetaI);
    Float thetaI = std::asin(clamp(sinThetaI, -1.0, 1.0));

    Float thetaD = (thetaI - thetaO) * 0.5;
    Float cosThetaD = cos(thetaD);

    Float thetaOR = thetaO - 2.0 * _scaleAngle;
    Float thetaOTT = thetaO + _scaleAngle;
    Float thetaOTRT = thetaO + 4 * _scaleAngle;

    Float MR = M(_vR, sin(thetaOR), sinThetaI, cos(thetaOR), cosThetaI);
    Float MTT = M(_vTT, sin(thetaOTT), sinThetaI, cos(thetaOTT), cosThetaI);
    Float MTRT = M(_vTRT, sin(thetaOTRT), sinThetaI, cos(thetaOTRT), cosThetaI);

    Float phi = std::atan2(event.wi.x, event.wi.z) - std::atan2(event.wo.x, event.wo.z);
    if (phi < 0.0f)
        phi += Constant::TWO_PI;

    Float h = getH(event);
    if (useGussianAmz) {
        return MR * _nR->eval(phi, cosThetaD)
               + MTT * _nTT->eval(phi, cosThetaD)
               + MTRT * _nTRT->eval(phi, cosThetaD);
    }
    vec3 Nr = vec3(NR(_betaR, trigInverse(event.wo.y), phi, h));
    vec3 Ntt = NP(_betaR, cosThetaD, phi, 1, h);
    vec3 Ntrt = NP(_betaR, cosThetaD, phi, 2, h);
    //  return Ntrt * MTRT;
    vec3 fsum = MR * Nr + MTT * Ntt + MTRT * Ntrt;
    //  return MR * Nr;
    return fsum;
    vec3 res = fsum;
    //  if( AbsCosTheta(event.wi)>0) res/= AbsCosTheta(event.wi);
    res = MTRT * Ntrt;
    return res;
}


Float Hair::Pdf(const SurfaceEvent &event) const {


    Float sinThetaO = event.wo.y;
    Float costhetaO = trigInverse(sinThetaO);
    Float thetaO = std::asin(clamp(sinThetaO, -1, 1));

    Float sinThetaI = event.wi.y;
    Float cosThetaI = sqrt(1 - sinThetaI * sinThetaI);
    Float thetaI = std::asin(clamp(sinThetaI, -1.0, 1.0));
    //First choose a lobe to sample
    Float h = getH(event);
    std::array<Float, pMax + 1> apPdf = ComputeApPdf(costhetaO, h);

    Float thetaOR = thetaO - 2.0 * _scaleAngle;
    Float thetaOTT = thetaO + _scaleAngle;
    Float thetaOTRT = thetaO + 4 * _scaleAngle;

    Float phiI = std::atan2(event.wi.x, event.wi.z);
    Float phiO = std::atan2(event.wo.x, event.wo.z);
    Float phi = phiI - phiO;

    Float cosThetaD = cos(0.5 * (std::asin(sinThetaI) - thetaO));
    Float iorPrime = std::sqrt(_eta * _eta - (1.0f - cosThetaD * cosThetaD)) / cosThetaD;
    Float gammaI = std::asin(clamp(h, -1.0f, 1.0f));
    Float gammaT = std::asin(clamp(h / iorPrime, -1.0f, 1.0f));

    Float pdf = 0;
    if (useGussianAmz) {
        phi = std::atan2(event.wi.x, event.wi.z);
        if (phi < 0.0f)
            phi += Constant::TWO_PI;

        float weightR = _nR->weight(cosThetaI);
        float weightTT = _nTT->weight(cosThetaI);
        float weightTRT = _nTRT->weight(cosThetaI);
        float weightSum = weightR + weightTT + weightTRT;

        float pdfR = weightR * M(_vR, sin(thetaOR), sinThetaI, cos(thetaOR), cosThetaI);
        float pdfTT = weightTT * M(_vTT, sin(thetaOTT), sinThetaI, cos(thetaOTT), cosThetaI);
        float pdfTRT = weightTRT * M(_vTRT, sin(thetaOTRT), sinThetaI, cos(thetaOTRT), cosThetaI);

        return (1.0f / weightSum) *
               (pdfR * _nR->pdf(phi, cosThetaD)
                + pdfTT * _nTT->pdf(phi, cosThetaD)
                + pdfTRT * _nTRT->pdf(phi, cosThetaD));
    }
    pdf += M(_vR, sin(thetaOR), sinThetaI, cos(thetaOR), cosThetaI) * apPdf[0] *
           TrimmedLogistic(phi - Phi(gammaI, gammaT, 0), _betaR, -Constant::PI, Constant::PI);
    pdf += M(_vTT, sin(thetaOTT), sinThetaI, cos(thetaOTT), cosThetaI) * apPdf[1] *
           TrimmedLogistic(phi - Phi(gammaI, gammaT, 1), _betaR, -Constant::PI, Constant::PI);
    pdf += M(_vTRT, sin(thetaOTRT), sinThetaI, cos(thetaOTRT), cosThetaI) * apPdf[2] *
           TrimmedLogistic(phi - Phi(gammaI, gammaT, 2), _betaR, -Constant::PI, Constant::PI);
    if (isnan(pdf)) {

    }
    return pdf;
}

Spectrum Hair::sampleF(SurfaceEvent &event, const vec2 &u) const {
    //Hair-samplineg requires 4 randoms.
    vec2 u0 = DemuxFloat(u[0]), u1 = DemuxFloat(u[1]);
    Float sinThetaO = event.wo.y;
    Float costhetaO = trigInverse(sinThetaO);
    Float thetaO = std::asin(clamp(sinThetaO, -1, 1));
    //First choose a lobe to sample
    Float h = getH(event);
    std::array<Float, pMax + 1> apPdf = ComputeApPdf(costhetaO, h);
    int p;
    for (p = 0; p < pMax; ++p) {
        if (u0[0] <= apPdf[p])
            break;
        u0[0] -= apPdf[p];
    }
    Float theta, v;
    if (p == 0) {
        theta = thetaO - 2 * _scaleAngle;
        v = _vR;
    } else if (p == 1) {
        theta = thetaO + _scaleAngle;
        v = _vTT;
    } else if (p == 2) {
        theta = thetaO + 4 * _scaleAngle;
        v = _vTRT;
    } else { throw ("Invalid P"); }

    Float sinThetaI = sampleM(v, sin(theta), cos(theta), u1[0], u1[1]);
    Float cosThetaI = trigInverse(sinThetaI);

    Float phiO = std::atan2(event.wo.x, event.wo.z);
    Float deltaphi;
    Float cosThetaD = cos(0.5 * (std::asin(clamp(sinThetaI, -1, 1)) - thetaO));
    Float iorPrime = std::sqrt(_eta * _eta - (1.0f - cosThetaD * cosThetaD)) / cosThetaD;
    Float gammaI = std::asin(clamp(h, -1.0f, 1.0f));
    Float gammaT = std::asin(clamp(h / iorPrime, -1.0f, 1.0f));

    Float phi;
    if (useGussianAmz) {
        Float phiPdf;
        if (p == 0)
            _nR->sample(cosThetaD, u0[1], phi, phiPdf);
        else if (p == 1)
            _nTT->sample(cosThetaD, u0[1], phi, phiPdf);
        else if (p == 2)
            _nTRT->sample(cosThetaD, u0[1], phi, phiPdf);
    } else {
        deltaphi = Phi(gammaI, gammaT, p) + SampleTrimmedLogistic(u0[1], _betaR, -Constant::PI, Constant::PI);
        phi = phiO + deltaphi;
    }

    Float sinPhi = sin(phi), cosPhi = cos(phi);

    event.wi = vec3(sinPhi * cosThetaI, sinThetaI, cosPhi * cosThetaI);
    event.sampleType = this->m_type;
    event.pdf = Pdf(event);
    return f(event);
}


Float Hair::D(Float beta, Float phi) const {
    Float result = 0.0f;
    Float delta;
    Float shift = 0.0f;
    do {
        delta = Gussian(beta, phi + shift) + Gussian(beta, phi - shift - Constant::TWO_PI);
        result += delta;
        shift += Constant::TWO_PI;
    } while (delta > 1e-4f);
    return result;
}

Hair::Hair(const Json &json) : BSDF(BXDFType(BSDF_GLOSSY | BSDF_TRANSMISSION | BSDF_REFLECTION)) {
    _scaleAngle = getOptional(json, "scale_angle", 2.5);
    Float melaninRatio = getOptional(json, "melanin_ratio", 1.f);
    Float melaninConcentration = getOptional(json, "melanin_concentration", 1.3);
    bool overrideSigmaA = containsAndGet(json, "sigma_a", _sigmaA);
    if (!overrideSigmaA) {
        const vec3 eumelaninSigmaA = vec3(0.419f, 0.697f, 1.37f);
        const vec3 pheomelaninSigmaA = vec3(0.187f, 0.4f, 1.05f);

        _sigmaA = melaninConcentration * lerp(eumelaninSigmaA, pheomelaninSigmaA, melaninRatio);
    }

    betaM = getOptional(json, "beta_m", 0.3);
    betaN = getOptional(json, "beta_n", 0.3);

    if (!useGussianAmz) {
        _betaR = std::max(Constant::PI * 0.5f * betaM, 0.04f);
        _betaTT = _betaR * 0.5f;
        _betaTRT = _betaR * 2.0f;

        _vR = _betaR * _betaR;
        _vTT = _betaTT * _betaTT;
        _vTRT = _betaTRT * _betaTRT;
    } else {
        _vR = sqr(0.726f * betaM + 0.812f * sqr(betaM) + 3.7f * std::pow(betaM, 20));
        _vTT = .25 * _vR;
        _vTRT = 4 * _vR;

        static const Float SqrtPiOver8 = 0.626657069f;
        _betaR = SqrtPiOver8 * (0.265f * betaN + 1.194f * sqr(betaN) + 5.372f * std::pow(betaN, 22));
    }


    _scaleAngle = Angle::degToRad(_scaleAngle);


    if (useGussianAmz)
        precomputeGussAmz();

}

//return sinTheta
Float Hair::sampleM(Float v, Float sinThetaO, Float cosThetaO, Float xi1, Float xi2) const {
    Float cosTheta = 1.0f + v * std::log(xi1 + (1.0f - xi1) * std::exp(-2.0f / v));
    Float sinTheta = trigInverse(cosTheta);
    Float cosPhi = std::cos(Constant::TWO_PI * xi2);
    return clamp(-cosTheta * sinThetaO + sinTheta * cosPhi * cosThetaO, -1, 1);
}

std::array<Float, pMax + 1> Hair::ComputeApPdf(Float cosThetaO, Float h) const {

    Float sinThetaO = trigInverse(cosThetaO);

    Float sinThetaT = sinThetaO / _eta;
    Float cosThetaT = trigInverse(sinThetaT);

    Float etap = std::sqrt(_eta * _eta - sqr(sinThetaO)) / cosThetaO;
    Float sinGammaT = h / etap;
    Float cosGammaT = trigInverse(sinGammaT);

    Spectrum T = exp(-_sigmaA * 2.f * cosGammaT / cosThetaT);
    std::array<Spectrum, pMax + 1> ap = Ap(cosThetaO, _eta, h, T);

    // Compute $A_p$ PDF from individual $A_p$ terms
    std::array<Float, pMax + 1> apPdf;
    Float sumY =
            std::accumulate(ap.begin(), ap.end() - 1, Float(0),
                            [](Float s, const Spectrum &ap) { return s + luminace(ap); });
    for (int i = 0; i < pMax; ++i) apPdf[i] = luminace(ap[i]) / sumY;
    if (isnan(apPdf[0])) {

    }
    return apPdf;
}

void Hair::precomputeGussAmz() {
    const int Resolution = PrecomputedAzimuthalLobe::AzimuthalResolution;
    std::unique_ptr<vec3[]> valuesR(new vec3[Resolution * Resolution]);
    std::unique_ptr<vec3[]> valuesTT(new vec3[Resolution * Resolution]);
    std::unique_ptr<vec3[]> valuesTRT(new vec3[Resolution * Resolution]);

    // Ideally we could simply make this a constexpr, but MSVC does not support that yet (boo!)
#define NumPoints 140

    GaussLegendre<NumPoints> integrator;
    const auto points = integrator.points();
    const auto weights = integrator.weights();

    // Cache the gammaI across all integration points
    std::array<float, NumPoints> gammaIs;
    for (int i = 0; i < NumPoints; ++i)
        gammaIs[i] = std::asin(points[i]);

    // Precompute the Gaussian detector and sample it into three 1D tables.
    // This is the only part of the precomputation that is actually approximate.
    // 2048 samples are enough to support the lowest roughness that the BCSDF
    // can reliably simulate
    const int NumGaussianSamples = 2048;
    std::unique_ptr<float[]> Ds[3];
    for (int p = 0; p < 3; ++p) {
        Ds[p].reset(new float[NumGaussianSamples]);
        for (int i = 0; i < NumGaussianSamples; ++i)
            Ds[p][i] = D(_betaR, i / (NumGaussianSamples - 1.0f) * Constant::TWO_PI);
    }

    // Simple wrapped linear interpolation of the precomputed table
    auto approxD = [&](int p, float phi) {
        float u = std::abs(phi * (Constant::INV_TWO_PI * (NumGaussianSamples - 1)));
        int x0 = int(u);
        int x1 = x0 + 1;
        u -= x0;
        return Ds[p][x0 % NumGaussianSamples] * (1.0f - u) + Ds[p][x1 % NumGaussianSamples] * u;
    };

    // Here follows the actual precomputation of the azimuthal scattering functions
    // The scattering functions are parametrized with the azimuthal angle phi,
    // and the cosine of the half angle, cos(thetaD).
    // This parametrization makes the azimuthal function relatively smooth and allows using
    // really low resolutions for the table (64x64 in this case) without any visual
    // deviation from ground truth, even at the lowest supported roughness setting
    for (int y = 0; y < Resolution; ++y) {
        float cosHalfAngle = y / (Resolution - 1.0f);

        // Precompute reflection Fresnel factor and reduced absorption coefficient
        float iorPrime = std::sqrt(_eta * _eta - (1.0f - cosHalfAngle * cosHalfAngle)) / cosHalfAngle;
        float cosThetaT = std::sqrt(1.0f - (1.0f - cosHalfAngle * cosHalfAngle) * sqr(1.0f / _eta));
        vec3 sigmaAPrime = _sigmaA / cosThetaT;

        // Precompute gammaT, f_t and internal absorption across all integration points
        std::array<float, NumPoints> fresnelTerms, gammaTs;
        std::array<vec3, NumPoints> absorptions;
        for (int i = 0; i < NumPoints; ++i) {
            gammaTs[i] = std::asin(clamp(points[i] / iorPrime, -1.0f, 1.0f));
            fresnelTerms[i] = Fresnel::dielectricReflectance(1.0f / _eta, cosHalfAngle * std::cos(gammaIs[i]));
            absorptions[i] = exp(-sigmaAPrime * 2.0f * std::cos(gammaTs[i]));
        }

        for (int phiI = 0; phiI < Resolution; ++phiI) {
            float phi = Constant::TWO_PI * phiI / (Resolution - 1.0f);

            float integralR = 0.0f;
            vec3 integralTT(0.0f);
            vec3 integralTRT(0.0f);

            // Here follows the integration across the fiber width, h.
            // Since we were able to precompute most of the factors that
            // are constant w.r.t phi for a given h,
            // we don't have to do much work here.
            for (int i = 0; i < integrator.numSamples(); ++i) {
                float fR = fresnelTerms[i];
                vec3 T = absorptions[i];

                float AR = fR;
                vec3 ATT = (1.0f - fR) * (1.0f - fR) * T;
                vec3 ATRT = ATT * fR * T;

                integralR += weights[i] * approxD(0, phi - Phi(gammaIs[i], gammaTs[i], 0)) * AR;
                integralTT += weights[i] * approxD(1, phi - Phi(gammaIs[i], gammaTs[i], 1)) * ATT;
                integralTRT += weights[i] * approxD(2, phi - Phi(gammaIs[i], gammaTs[i], 2)) * ATRT;
            }

            valuesR[phiI + y * Resolution] = vec3(0.5f * integralR);
            valuesTT[phiI + y * Resolution] = 0.5f * integralTT;
            valuesTRT[phiI + y * Resolution] = 0.5f * integralTRT;
        }
    }

    // Hand the values off to the helper class to construct sampling CDFs and so forth
    _nR.reset(new PrecomputedAzimuthalLobe(std::move(valuesR)));
    _nTT.reset(new PrecomputedAzimuthalLobe(std::move(valuesTT)));
    _nTRT.reset(new PrecomputedAzimuthalLobe(std::move(valuesTRT)));
}




