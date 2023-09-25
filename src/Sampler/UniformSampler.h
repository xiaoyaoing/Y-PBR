#include "Sampler.hpp"
#include "Rng.h"
class UniformSampler : public Sampler {
public:
    UniformSampler(int spp = 4,uint64_t seed =  0) ;

    std::unique_ptr<Sampler> clone(int seed) const override ;

    Float getNext1D() override { return rng.getNext(); }
    vec2 getNext2D() override { return vec2 (rng.getNext(), rng.getNext()); }

protected:
    RNG rng;
};