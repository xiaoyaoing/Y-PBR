#include "Sampler.hpp"
#include "Rng.h"
class UniformSampler : public Sampler {
public:
    UniformSampler(int spp = 4,int seed =  0) : Sampler(spp),rng(seed) {}

    std::unique_ptr<Sampler> clone(int seed) const override {
        return std::make_unique <UniformSampler>(samplesPerPixel,seed);
    }

    Float getNext1D() override { return rng.getNext(); }
    vec2 getNext2D() override { return vec2 (rng.getNext(), rng.getNext()); }

protected:
    RNG rng;
};