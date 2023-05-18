#include "Render.hpp"

#include "Common/Debug.hpp"
#include "IO/FileUtils.hpp"
#include "Debug.hpp"

#include <Integrator/PhotonMapper.hpp>
#include <Integrator/VolPathIntegrator.h>
#include <Integrator/BDPT/BDPTIntegrator.hpp>
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
        const std::string & type = getOptional(integratorJson, "type", std::string("path_tracer"));
        if ( type == "path_tracer" ) {
            std::string lightSampleStrategy = getOptional(integratorJson, "light_sample", std::string("uniform"));
            integrator = std::make_unique < PathIntegrator >(camera, sampler, integratorJson);
        } else if ( type == "vol_path_tracer" ) {
            std::string lightSampleStrategy = getOptional(integratorJson, "light_sample", std::string("uniform"));
            integrator = std::make_unique < VolPathIntegrator >(camera, sampler, integratorJson);

        } else if ( type == "sppm" ) {
            int maxBounces = getOptional(integratorJson, "max_bounces", 8);
            integrator = std::make_unique < PhotonMapper >(camera, integratorJson);
        }
        else if(type == "bidirectional_path_tracer"){
            integrator = std::make_unique < BDPTIntegrator >(camera, sampler,integratorJson);
        }
    }
}

class TimeCounter{
public:
    void printTimeCount(){
        std::chrono::time_point < std::chrono::steady_clock > end = std::chrono::steady_clock::now();
        auto s = std::chrono::duration_cast < std::chrono::seconds >(end - start).count();
        std::cout <<std::endl<< event<<"done. Take " << s << "s";
    }
    TimeCounter(std::string event):event(std::move(event))
    {
        start = std::chrono::steady_clock::now();
    }
private:
    std::chrono::time_point < std::chrono::steady_clock >  start;
    std::string event;

};

void Render::Go( ) {
    TimeCounter counter("Rendering");
    scene->build();
    integrator->process(* scene, * sampler);
    integrator->render(* scene);
    counter.printTimeCount();
}