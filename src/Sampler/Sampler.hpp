#ifndef _SAMPLER_H
#define _SAMPLER_H
#include <cstdint>
#include <limits>
#include <memory>

#include "Common/math.hpp"

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
typedef struct {
    uint64_t state;
    uint64_t inc;
} pcg32_random_t;

inline uint32_t pcg32_random_r(pcg32_random_t* rng) {
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// random number generator
class RNG {
private:
    pcg32_random_t state;

public:
    RNG() {
        state.state = 1;
        state.inc = 1;
    }
    RNG(uint64_t seed) {
        state.state = seed;
        state.inc = 1;
    }

    uint64_t getSeed() const { return state.state; }
    void setSeed(uint64_t seed) { state.state = seed; }

    Float getNext() {
        constexpr Float divider = 1.0 / std::numeric_limits<uint32_t>::max();
        return pcg32_random_r(&state) * divider;
    }
};

// sampler interface
class Sampler {
protected:
    RNG rng;

public:
    Sampler() {}
    virtual  ~Sampler() = default;
    Sampler(uint64_t seed) : rng(seed) {}

    uint64_t getSeed() const { return rng.getSeed(); }
    void setSeed(uint64_t seed) { rng.setSeed(seed); }

    virtual std::unique_ptr<Sampler> clone() const = 0;
    virtual Float getNext1D() = 0;
    virtual vec2 getNext2D() = 0;
};

// uniform distribution sampler
class UniformSampler : public Sampler {
public:
    UniformSampler() : Sampler() {}
    UniformSampler(uint64_t seed) : Sampler(seed) {}

    std::unique_ptr<Sampler> clone() const override {
        return std::make_unique <UniformSampler>();
    }

    Float getNext1D() override { return rng.getNext(); }
    vec2 getNext2D() override { return vec2 (rng.getNext(), rng.getNext()); }
};

class HaltonSampler : public Sampler{

};

//// sample direction in the hemisphere
//// its pdf is propotional to cosine
//inline Vec3f sampleCosineHemisphere(const Vec2f& uv, float& pdf) {
//    const float theta =
//            0.5f * std::acos(std::clamp(1.0f - 2.0f * uv[0], -1.0f, 1.0f));
//    const float phi = PI_MUL_2 * uv[1];
//    const float cosTheta = std::cos(theta);
//    pdf = PI_INV * cosTheta;
//    return sphericalToCartesian(theta, phi);
//}


//struct Distribution1D {
//    // Distribution1D Public Methods
//    Distribution1D(const Float *f, int n) : func(f, f + n), cdf(n + 1) {
//        // Compute integral of step function at $x_i$
//        cdf[0] = 0;
//        for (int i = 1; i < n + 1; ++i) cdf[i] = cdf[i - 1] + func[i - 1] / n;
//
//        // Transform step function integral into CDF
//        funcInt = cdf[n];
//        if (funcInt == 0) {
//            for (int i = 1; i < n + 1; ++i) cdf[i] = Float(i) / Float(n);
//        } else {
//            for (int i = 1; i < n + 1; ++i) cdf[i] /= funcInt;
//        }
//    }
//    int Count() const { return (int)func.size(); }
//    Float SampleContinuous(Float u, Float *pdf, int *off = nullptr) const {
//        // Find surrounding CDF segments and _offset_
//        int offset = FindInterval((int)cdf.size(),
//                                  [&](int index) { return cdf[index] <= u; });
//        if (off) *off = offset;
//        // Compute offset along CDF segment
//        Float du = u - cdf[offset];
//        if ((cdf[offset + 1] - cdf[offset]) > 0) {
//            CHECK_GT(cdf[offset + 1], cdf[offset]);
//            du /= (cdf[offset + 1] - cdf[offset]);
//        }
//        DCHECK(!std::isnan(du));
//
//        // Compute PDF for sampled offset
//        if (pdf) *pdf = (funcInt > 0) ? func[offset] / funcInt : 0;
//
//        // Return $x\in{}[0,1)$ corresponding to sample
//        return (offset + du) / Count();
//    }
//    int SampleDiscrete(Float u, Float *pdf = nullptr,
//                       Float *uRemapped = nullptr) const {
//        // Find surrounding CDF segments and _offset_
//        int offset = FindInterval((int)cdf.size(),
//                                  [&](int index) { return cdf[index] <= u; });
//        if (pdf) *pdf = (funcInt > 0) ? func[offset] / (funcInt * Count()) : 0;
//        if (uRemapped)
//            *uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
//        if (uRemapped) CHECK(*uRemapped >= 0.f && *uRemapped <= 1.f);
//        return offset;
//    }
//    Float DiscretePDF(int index) const {
//        CHECK(index >= 0 && index < Count());
//        return func[index] / (funcInt * Count());
//    }
//
//    // Distribution1D Public Data
//    std::vector<Float> func, cdf;
//    Float funcInt;
//};

#endif