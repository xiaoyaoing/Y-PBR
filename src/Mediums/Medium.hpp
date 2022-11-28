#include "Ray/Ray.hpp"
#include "Colors/Spectrum.hpp"
#include "Sampler/Sampler.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"

class  Medium {
public:
    virtual Float TR(const Ray & ray) const =  0;
    virtual Spectrum sampleDistance(const Ray & ray,Sampler & sampler,VolumeEvent & event) const = 0;
};

class Homogeneous : public  Medium {
    Float sigmaT = 1.5;
public:
    Float TR(const Ray & ray) const override;

    Spectrum sampleDistance(const Ray & ray, Sampler & sampler, VolumeEvent & event) const override;
};
