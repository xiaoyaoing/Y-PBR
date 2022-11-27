#include "VolPathIntegrator.h"

void VolPathIntegrator::process(const Scene & scene, Sampler & sampler) {

}

vec3 VolPathIntegrator::integrate(const Ray & ray, const Scene & scene, Sampler & sampler) const {
    const Medium * medium;
    medium = _camera->Medium();
    return vec3();
}
