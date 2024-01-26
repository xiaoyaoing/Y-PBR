#pragma once
#include <cstdint>
#include "limits"
// random number generator
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
    uint32_t rot        = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

class RNG {
private:
    pcg32_random_t         state;
    static constexpr float divider = 1.0 / std::numeric_limits<uint32_t>::max();

public:
    RNG() {
        state.state = 1;
        state.inc   = 1;
    }
    RNG(uint64_t seed) {
        state.state = seed;
        state.inc   = 1;
    }

    uint64_t getSeed() const { return state.state; }
    void     setSeed(uint64_t seed) { state.state = seed; }

    float getNext() {
        return pcg32_random_r(&state) * divider;
    }

    uint32_t UniformUInt32(uint32_t b) {
        uint32_t threshold = (~b + 1u) % b;
        while (true) {
            uint32_t r = pcg32_random_r(&state);
            if (r >= threshold) return r % b;
        }
    }
};