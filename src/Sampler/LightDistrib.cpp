#include "LightDistrib.hpp"

//todo power distribution
std::unique_ptr<Distribution1D> CreateLightSampleDistribution(
    const std::string& name,
    const Scene&       scene) {
    std::vector<Float> lightWeights;
    int                nLights = scene.lights.size();
    if (name == "uniform") {
        lightWeights.resize(nLights);
        std::fill(lightWeights.begin(), lightWeights.end(), 1);
    } else if (name == "power") {
        for (const auto& light : scene.lights)
            lightWeights.push_back(luminace(light->Power()));
    }
    return std::make_unique<Distribution1D>(lightWeights.data(), nLights);
}

UniformLightDistribution::UniformLightDistribution(const Scene& scene) {
    std::vector<Float> prob(scene.lights.size(), Float(1));
    distrib.reset(new Distribution1D(&prob[0], int(prob.size())));
}
const Distribution1D* UniformLightDistribution::Lookup(const vec3& p) const {
    return distrib.get();
}