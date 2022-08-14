#pragma once
#include "Distrib.hpp"
#include "../scene.hpp"

class LightDistribution {
public:
    virtual ~LightDistribution() = default;;

    // Given a point |p| in space, this method returns a (hopefully
    // effective) sampling distribution for light sources at that point.
    virtual const Distribution1D *Lookup(const vec3 &p) const = 0;
};

class UniformLightDistribution : public  LightDistribution{
public:
    UniformLightDistribution(const Scene &scene);
    const Distribution1D *Lookup(const vec3 &p) const;

private:
    std::unique_ptr<Distribution1D> distrib;
};

std::unique_ptr<LightDistribution> CreateLightSampleDistribution(
        const std::string &name, const Scene &scene);



