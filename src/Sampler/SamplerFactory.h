#include "Sampler.hpp"
#include "Common/Json.hpp"
namespace SamplerFactory {
    std::shared_ptr<Sampler> loadSampler(const Json& json, int spp, ivec2 res);
}