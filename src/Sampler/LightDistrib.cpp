#include "LightDistrib.hpp"


//todo power distribution
std::unique_ptr<LightDistribution> CreateLightSampleDistribution(
        const std::string &name, const Scene &scene){
    return std::make_unique<UniformLightDistribution>(scene);
    if(name=="uniform"){
        return std::make_unique<UniformLightDistribution>(scene);
    }
}

UniformLightDistribution::UniformLightDistribution(const Scene &scene){
    std::vector<Float> prob(scene.lights.size(), Float(1));
    distrib.reset(new Distribution1D(&prob[0], int(prob.size())));
}
const Distribution1D * UniformLightDistribution::Lookup(const vec3 &p) const{
    return distrib.get();
}






