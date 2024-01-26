#pragma once

#include "Sampler.hpp"
#include "LowDiscrepancy.hpp"
#include <vector>

// HaltonSampler Declarations
class HaltonSampler : public GlobalSampler {
public:
    // HaltonSampler Public Methods
    HaltonSampler(int nsamp, ivec4 sampleBounds, bool sampleAtCenter = false);

    int64_t GetIndexForSample(int64_t sampleNum) const;

    Float SampleDimension(int64_t index, int dimension) const;

    std::unique_ptr<Sampler> Clone(int seed);

    std::unique_ptr<Sampler> clone(int seed) const override;

private:
    // HaltonSampler Private Data
    static std::vector<uint16_t> radicalInversePermutations;
    ivec2                        baseScales, baseExponents;
    int                          sampleStride;
    int                          multInverse[2];
    mutable ivec2                pixelForOffset = ivec2(std::numeric_limits<int>::max(),
                                         std::numeric_limits<int>::max());
    mutable int64_t              offsetForCurrentPixel;
    // Added after book publication: force all image samples to be at the
    // center of the pixel area.

    // HaltonSampler Private Methods
    const uint16_t* PermutationForDimension(int dim) const {
        if (dim >= PrimeTableSize)
            std::runtime_error("Dim excedd limit");
        return &radicalInversePermutations[PrimeSums[dim]];
    }
};