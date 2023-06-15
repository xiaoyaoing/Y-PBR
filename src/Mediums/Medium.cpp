#include "Medium.hpp"
#include "PhaseFunctionFactory.hpp"

bool GetMediumScatteringProperties(const std::string &name, Spectrum *sigma_a,
                                   Spectrum *sigma_prime_s) {
    for (MeasuredSS &mss : SubsurfaceParameterTable) {
        if (name == mss.name) {
            *sigma_a = mss.sigma_a;
            *sigma_prime_s = mss.sigma_prime_s;
            return true;
        }
    }
    return false;
}

Spectrum Homogeneous::TR(const Ray & ray) const {
    return  exp((ray.nearT - ray.farT) * sigmaT);
}

Spectrum Homogeneous::sampleDistance(const Ray & ray, Sampler & sampler, VolumeEvent & event) const {
    if(isBlack(sigmaS)){
        auto maxT = ray.farT;
        if (maxT == std::numeric_limits<Float>::max())
            return Spectrum (0);
        event.p = ray.operator()(maxT);
        event.phase = phaseFunction.get();
        event.exited = true;
        event.rayDir = ray.d;
        event.pdf = 1.0f;
        return exp(-maxT * sigmaT);
    }
    int channel = std::min(int(sampler.getNext1D() * 3),2);
    Float t = -std::log(1-sampler.getNext1D()) /sigmaT[channel];
    bool hitSurace = t >= ray.farT;
    t = std::min(t,ray.farT);
    event.exited =hitSurace;
    event.phase = phaseFunction.get();
    event.p = ray(t);
    event.rayDir = ray.d;
    Spectrum Tr = exp(-t * sigmaT);
    event.pdf = average(hitSurace?Tr:Tr * sigmaT);
   // return Tr;
    return  (hitSurace?Tr/event.pdf:Tr * sigmaS / event.pdf);
}

Homogeneous::Homogeneous(const Json & json) {
    phaseFunction = PhaseFunctionFactory::loadPhaseFromJson(json.at("phase_function"));
    Float scale = getOptional<Float>(json,"density",Float(1));
    sigmaA = getOptional(json,"sigma_a",Spectrum(0.5)) * scale;
    sigmaS = getOptional(json,"sigma_s",Spectrum(0.5)) * scale;
    sigmaT = sigmaA + sigmaS;
}
