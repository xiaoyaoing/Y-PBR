#include "SamplerFactory.h"
#include "HaltonSampler.h"
#include "Sampler/UniformSampler.h"

std::shared_ptr<Sampler> SamplerFactory::loadSampler(const Json &json,int spp,ivec2 res) {
   // std::string type = getOptional(json,"type","uniform");
   std::string type = getOptional(json,"type",std::string("uniform"));
   if(type == "uniform")
       return std::make_shared<UniformSampler>(spp);
   if(type == "halton")
       return std::make_shared<HaltonSampler>(spp,ivec4(0,0,res));
    return std::shared_ptr<Sampler>();
}
