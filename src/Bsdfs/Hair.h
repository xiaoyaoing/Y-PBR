#include "Reflection.hpp"
#include "MicrofacetDistribution.hpp"
#include "PrecomputeALobe.h"

const int pMax = 3;

/// Hair Bsdf
class Deielectric;
class Hair : public  BSDF{
public:
    Spectrum f(const SurfaceEvent & event) const override;

    Float Pdf(const SurfaceEvent & event) const override;

    Spectrum sampleF(SurfaceEvent & event, const vec2 & u) const override;

    Hair(const Json & json);
protected:
    Float M(Float v,Float sinThetaO,Float sinThetaI,Float cosThetaO,Float cosThetaI) const ;
    Float sampleM(Float v,Float sinThetaO,Float cosThetaO,Float xi1,Float xi2 ) const ;

    Float NR(Float beta,Float cosO,Float phi,Float h) const;
    vec3 NP(Float beta,Float cosO,Float phi,int p,Float h) const;
    Float D(Float beta,Float phi) const;

    std::array<Float, pMax + 1> ComputeApPdf(Float cosThetaO,Float h) const;

    void precomputeGussAmz();

private:
    std::unique_ptr<PrecomputedAzimuthalLobe> _nR = nullptr, _nTT = nullptr, _nTRT = nullptr;

    vec2 alpha_r,alpha_tt,alpha_trt;
    Float betaM,betaN;
    Float _scaleAngle;
    Float _roughness;
    Float _betaR,_betaTT,_betaTRT;
    Float _vR,_vTT,_vTRT;
    vec3 _sigmaA;
    const Float _eta = 1.55;
    bool useGussianAmz ;
};

// https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
static uint32_t Compact1By1(uint32_t x) {
    // TODO: as of Haswell, the PEXT instruction could do all this in a
    // single instruction.
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x &= 0x55555555;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >> 1)) & 0x33333333;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >> 2)) & 0x0f0f0f0f;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >> 4)) & 0x00ff00ff;
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x >> 8)) & 0x0000ffff;
    return x;
}

static vec2 DemuxFloat(Float f) {
    uint64_t v = f * (1ull << 32);
    uint32_t bits[2] = {Compact1By1(v), Compact1By1(v >> 1)};
    return {bits[0] / Float(1 << 16), bits[1] / Float(1 << 16)};
}