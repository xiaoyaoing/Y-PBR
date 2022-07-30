#include "PathIntegrator.hpp"
//2022/7/15
vec3 PathIntegrator::integrate(const Ray &ray, const Scene &scene, Sampler &sampler) const {
    Intersection intersection;
    if(!scene.intersect(ray,intersection)){
        return vec3(0);
    }

    return  intersection.normal;
    vec3 lightPos(1,1,1);
    vec3 lightPower(0.1,0.5,0.5);
    vec3 dir = lightPos - intersection.pos;
    return lightPower / length2(dir) * dot(normalize(dir),intersection.normal);

}

PathIntegrator::PathIntegrator(nlohmann::json j) : Integrator(j) {

}
