#include "UniformSampler.h"

UniformSampler::UniformSampler(int spp, uint64_t seed) : Sampler(spp)
{
    rng.setSeed(seed);
}

std::unique_ptr<Sampler> UniformSampler::clone(int seed) const {

        return std::make_unique <UniformSampler>(samplesPerPixel,seed);
    }
