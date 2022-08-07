#include "PathIntegrator.hpp"

//2022/7/15
Spectrum PathIntegrator::integrate(const Ray &ray, const Scene &scene, Sampler &sampler) const {
//    Intersection intersection =


//    vec3 lightPos(-20, 40, 20);
//    vec3 lightPower(2e3, 2.e3, 2e3);
//    vec3 dir = lightPos - intersection.pos;
//    return lightPower / length2(dir);

     std::optional<Intersection> its ;
     vec3 val;
     Spectrum throughPut(1.0);
     Spectrum  L;
     int bounces=0,maxDepth=10;

    for (bounces = 0;; ++bounces) {
         its = scene.intersect(ray);
         if (!its.has_value() || bounces >= maxDepth) break;

        // return  glm::abs(its->n) ;

         L += throughPut * its->Le(-ray.d);
         its->shFrame=Frame(its->n);
         its->wo= normalize(its->shFrame.toLocal(-ray.d));
         auto  Ld = UniformSampleOneLight(its.value(),scene,sampler);  //direct lighting


         L+=Ld;
         return Ld;

//         auto wo = its.shFrame.toLocal(-ray.d);
//
//         BsdfSample bsdfSample;
//         its->bsdf->sample(bsdfSample);
//
//         throughPut *= bsdfSample.val * cositem;
//         ray = Ray(its->pos,dir);
     }
    return Spectrum(0);



}

PathIntegrator::PathIntegrator(nlohmann::json j) : Integrator(j) {

}
