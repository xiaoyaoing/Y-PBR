#include "Medium.hpp"
#include "PhaseFunctionFactory.hpp"
Spectrum Homogeneous::TR(const Ray & ray) const {
    return  exp((ray.nearT - ray.farT) * sigmaT);
}

Spectrum Homogeneous::sampleDistance(const Ray & ray, Sampler & sampler, VolumeEvent & event) const {
    int channel = std::min(int(sampler.getNext1D() * 3),2);
    Float t = -std::log(1-sampler.getNext1D()) /sigmaT[channel];
    bool hitSurace = t >= ray.farT;
    t = std::min(t,ray.farT);
    event.exited =hitSurace;
    event.phase = phaseFunction.get();
    event.p = ray(t);
    event.rayDir = ray.d;
    Spectrum Tr = exp(-t * sigmaT);
    event.pdf = average(hitSurace?Tr:Tr/sigmaT);
    //return Spectrum(1);
    return  hitSurace?Tr/event.pdf:Tr * sigmaS / event.pdf;
}

Homogeneous::Homogeneous(const Json & json) {
    phaseFunction = PhaseFunctionFactory::loadPhaseFromJson(json.at("phase_function"));
    Float scale = getOptional<Float>(json,"density",Float(1));
    sigmaA = getOptional(json,"sigma_a",Spectrum(0.5)) * scale;
    sigmaS = getOptional(json,"sigma_s",Spectrum(0.5)) * scale;
    sigmaT = sigmaA + sigmaS;
}
