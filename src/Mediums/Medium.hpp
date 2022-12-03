#pragma once
#include "Ray/Ray.hpp"
#include "Colors/Spectrum.hpp"
#include "Sampler/Sampler.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"

class  Medium {
public:
    virtual Spectrum TR(const Ray & ray) const =  0;
    virtual Spectrum sampleDistance(const Ray & ray,Sampler & sampler,VolumeEvent & event) const = 0;
};

class Homogeneous : public  Medium {
public:
    Homogeneous(const Json & json);
    Spectrum TR(const Ray & ray) const override;
    Spectrum sampleDistance(const Ray & ray, Sampler & sampler, VolumeEvent & event) const override;
protected:
    std::shared_ptr<PhaseFunction> phaseFunction;
    Spectrum sigmaA, sigmaS, sigmaT;

};
