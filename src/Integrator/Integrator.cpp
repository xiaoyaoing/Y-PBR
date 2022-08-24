#include "Integrator.hpp"
#include "iostream"
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





Spectrum
Integrator::EstimateDirect(const Intersection & its, const vec2 & uShading,
                           const Light & light, const vec2 & uLight,
                           const Scene & scene, Sampler & sampler, bool specular) const{

    BXDFType bsdfFlags =
            specular ? BSDF_ALL : BXDFType(BSDF_ALL & ~BSDF_SPECULAR);
    vec3 wi;
    Float lightPdf = 0, scatteringPdf = 0;
    VisibilityTester visibility;

    Spectrum  Ld(0);

    auto Li = light.Sample_Li(its,uLight,&wi,&lightPdf,&visibility);
//
    if(isBlack(Li)
    // || !visibility.Unoccluded(scene)
    ){
        return Spectrum();
    }

    if(lightPdf>0)
    {   auto f = its.bsdf->f(its.wo, its.toLocal(-wi),bsdfFlags) * abs(dot(wi, its.getNormal()));
        auto t=its.toLocal(-wi);
       // f=vec3 (1,1,1);
        Ld+= Li * f / lightPdf;
        if( isBlack(Ld) && its.bsdf->name=="shortBox"){
            std::string s = its.bsdf->name;
            int k=1;
        }
    }

    return Ld;
}

Spectrum
Integrator::UniformSampleOneLight(const Intersection & its,
                                  const Scene & scene,
                                  Sampler & sampler,
                                  const Distribution1D *lightDistrib,
                                  bool handleMedia) const {

    //todo multiple lightPdf
    Float  lightPdf;
    std::shared_ptr<Light> light;
    int lightNum;
    if (lightDistrib) {
        lightNum = lightDistrib->SampleDiscrete(sampler.getNext1D(), &lightPdf);
        if (lightPdf == 0) return Spectrum(0.f);
    } else {
        int nLights = scene.lights.size();
        lightNum = std::min((int)(sampler.getNext1D() * nLights), nLights - 1);
        lightPdf = Float(1) / nLights;
    }
    light = scene.lights.at(lightNum);

    vec2 uScattering = sampler.getNext2D();
    vec2 uLight=sampler.getNext2D();
    return EstimateDirect(its, uScattering, *light, uLight,
                          scene, sampler,  handleMedia) / lightPdf;
}
