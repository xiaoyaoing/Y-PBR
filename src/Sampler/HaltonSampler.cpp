#include "HaltonSampler.h"
static int kMaxResolution = 128;
static void extendedGCD(uint64_t a, uint64_t b, int64_t *x, int64_t *y) {
    if (b == 0) {
        *x = 1;
        *y = 0;
        return;
    }
    int64_t d = a / b, xp, yp;
    extendedGCD(b, a % b, &xp, &yp);
    *x = yp;
    *y = xp - (d * yp);
}
static uint64_t multiplicativeInverse(int64_t a, int64_t n) {
    int64_t x, y;
    extendedGCD(a, n, &x, &y);
    return mod(x, n);
}
std::vector<uint16_t > HaltonSampler::radicalInversePermutations;
HaltonSampler::HaltonSampler(int nsamp,  ivec4 sampleBounds, bool sampleAtCenter) : GlobalSampler(nsamp){
    if (radicalInversePermutations.empty()) {
        RNG rng;
        radicalInversePermutations = ComputeRadicalInversePermutations(rng);
    }
    // Find radical inverse base scales and exponents that cover sampling area
    ivec2 res = ivec2(sampleBounds.z,sampleBounds.w) -  ivec2(sampleBounds.x,sampleBounds.y);
    for (int i = 0; i < 2; ++i) {
        int base = (i == 0) ? 2 : 3;
        int scale = 1, exp = 0;
        while (scale < std::min(res[i], kMaxResolution)) {
            scale *= base;
            ++exp;
        }
        baseScales[i] = scale;
        baseExponents[i] = exp;
    }

    // Compute stride in samples for visiting each pixel area
    sampleStride = baseScales[0] * baseScales[1];

    // Compute multiplicative inverses for _baseScales_
    multInverse[0] = multiplicativeInverse(baseScales[1], baseScales[0]);
    multInverse[1] = multiplicativeInverse(baseScales[0], baseScales[1]);
}

int64_t HaltonSampler::GetIndexForSample(int64_t sampleNum) const {
    if (curPixel != pixelForOffset) {
        // Compute Halton sample offset for _currentPixel_
        offsetForCurrentPixel = 0;
        if (sampleStride > 1) {
            ivec2 pm(mod(curPixel[0], kMaxResolution),
                       mod(curPixel[1], kMaxResolution));
            for (int i = 0; i < 2; ++i) {
                uint64_t dimOffset =
                        (i == 0)
                        ? InverseRadicalInverse<2>(pm[i], baseExponents[i])
                        : InverseRadicalInverse<3>(pm[i], baseExponents[i]);
                offsetForCurrentPixel +=
                        dimOffset * (sampleStride / baseScales[i]) * multInverse[i];
            }
            offsetForCurrentPixel %= sampleStride;
        }
        pixelForOffset = curPixel;
    }
    return offsetForCurrentPixel + sampleNum * sampleStride;}

Float HaltonSampler::SampleDimension(int64_t index, int dimension) const {
   // if (sampleAtPixelCenter && (dim == 0 || dim == 1)) return 0.5f;
    if (dimension == 0)
        return RadicalInverse(dimension, index >> baseExponents[0]);
    else if (dimension == 1)
        return RadicalInverse(dimension, index / baseScales[1]);
    else
        return ScrambledRadicalInverse(dimension, index,
                                       PermutationForDimension(dimension));}

std::unique_ptr<Sampler> HaltonSampler::Clone(int seed) {
    return std::make_unique<HaltonSampler>(*this);
    //return std::make_unique<Sampler>(*this);
}

std::unique_ptr<Sampler> HaltonSampler::clone(int seed) const {
    return std::make_unique<HaltonSampler>(*this);
}
