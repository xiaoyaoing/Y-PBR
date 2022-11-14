#include "Render.hpp"

#include "Common/Debug.hpp"
#include "IO/FileUtils.hpp"
#include "Debug.hpp"

#include <Integrator/PhotonMapper.hpp>
#include <Integrator/VolPathIntegrator.h>
#include <thread>
#include <spdlog/spdlog.h>

#include <iostream>


Render::Render(const Json & json) {
    camera = std::make_unique < Camera >(json.at("camera"));
    scene = std::make_unique < Scene >(json);
    sampler = std::make_shared < UniformSampler >();

    //load integrator
    {
       const Json & integratorJson = json.at("integrator");
       const std::string & type = getOptional(integratorJson,"type",std::string("path_tracer"));
       if(type == "path_tracer")
       {
           std::string lightSampleStrategy = getOptional(integratorJson,"light_sample",std::string("uniform"));
           integrator = std::make_unique<PathIntegrator>(camera,sampler,lightSampleStrategy);
       }
       else if(type == "sppm")
       {
            Float radius = getOptional(integratorJson,"radius",0.05);
            int iterations = getOptional(integratorJson,"interation_num",64);
            int photonsPerIteration = getOptional(integratorJson,"photons_per",camera->image->width() * camera->image->height());
            int writeFrequency = getOptional(integratorJson,"write_frequency",32);
            int maxBounces = getOptional(integratorJson,"max_bounces",8);
            integrator = std::make_unique<PhotonMapper>(camera,iterations,radius,maxBounces,photonsPerIteration,writeFrequency);
       }
    }
}


void Render::Go( ) {
    std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
    scene->build();
    integrator->process(* scene, * sampler);
    integrator->render(*scene);

    std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();
    auto s =std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
    std::cout<<"\nrendering done take "<<s<<"s";
}