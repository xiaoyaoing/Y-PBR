#include "Integrator.hpp"

Integrator::Integrator(nlohmann::json j)
{

}


//Spectrum Integrator::sampleDirect(const Intersection& intersection, const Scene & scene,LightSample& ls) const
//{
//    if (scene.lights.empty()){
//
//    }
//
//
//    auto u = Sampler::get<Dim::LIGHT, 3>();
//
//    // Pick one light source and divide with probability of selecting light source
//    ls.light = scene.selectLight(u[2], ls.select_probability);
//
//    glm::dvec3 light_pos = ls.light->operator()(u[0], u[1]);
//    Ray shadow_ray(interaction.position + interaction.normal * C::EPSILON, light_pos);
//
//    double cos_light_theta = glm::dot(-shadow_ray.direction, ls.light->normal(light_pos));
//
//    if (cos_light_theta <= 0.0)
//    {
//        return glm::dvec3(0.0);
//    }
//
//    double cos_theta = glm::dot(shadow_ray.direction, interaction.normal);
//    if (cos_theta <= 0.0)
//    {
//        if (interaction.material->opaque || cos_theta == 0.0)
//        {
//            return glm::dvec3(0.0);
//        }
//        else
//        {
//            // Try transmission
//            shadow_ray = Ray(interaction.position - interaction.normal * C::EPSILON, light_pos);
//        }
//    }
//
//    Intersection shadow_intersection = scene.intersect(shadow_ray);
//
//    if (!shadow_intersection || shadow_intersection.surface != ls.light)
//    {
//        return glm::dvec3(0.0);
//    }
//
//    double light_pdf = pow2(shadow_intersection.t) / (ls.light->area() * cos_light_theta);
//
//    double bsdf_pdf;
//    glm::dvec3 bsdf_absIdotN;
//    if (!interaction.BSDF(bsdf_absIdotN, shadow_ray.direction, bsdf_pdf))
//    {
//        return glm::dvec3(0.0);
//    }
//
//    double mis_weight = powerHeuristic(light_pdf, bsdf_pdf);
//
//    return mis_weight * bsdf_absIdotN * ls.light->material->emittance / (light_pdf * ls.select_probability);
//}

vec3 Integrator::sampleDirectLight(const Scene & scene, const Intersection & intersection, Light sample) {

//    //todo mult light
//
//    auto light = scene.lights[0];
//
//    auto Le=light->sampleLi(intersection, sample,
//                vec3 * wi, Float *pdf,
//            VisibilityTester *vis);
//
//    return light;
}



Spectrum
Integrator::EstimateDirect(const Intersection & its, const vec2 & uShading, const Light & light, const vec2 & uLight,
                           const Scene & scene, Sampler & sampler, bool specular) const{


    vec3 wi;
    Float lightPdf = 0, scatteringPdf = 0;
    VisibilityTester visibility;
    Spectrum  Ld(0.0);

    auto Li = light.Sample_Li(its,uLight,&wi,&lightPdf,&visibility);

    auto f = its.bsdf->f(its.wo, its.shFrame.toLocal(wi)) * abs(dot(wi, its.n));

    return Li * f / lightPdf;

}

Spectrum
Integrator::UniformSampleOneLight(const Intersection & its, const Scene & scene, Sampler & sampler, bool handleMedia) const {

    //todo multiple lightPdf
    Float  lightPdf =1.0;
    auto light = scene.lights[0];

    vec2 uScattering = sampler.getNext2D();
    vec2 uLight=sampler.getNext2D();

    return EstimateDirect(its, uScattering, *light, uLight,
                          scene, sampler,  handleMedia) / lightPdf;
}
