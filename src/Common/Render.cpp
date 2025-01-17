#include "Render.hpp"

#include "Common/Debug.hpp"
#include "IO/FileUtils.hpp"
#include "Debug.hpp"

#include <Integrator/PhotonMapper.hpp>
#include <Integrator/VolPathIntegrator.h>
#include <Integrator/BDPT/BDPTIntegrator.hpp>
#include <thread>

#include <Sampler/UniformSampler.h>

#include <iostream>
#include "Sampler/SamplerFactory.h"
Render::Render(const Json& json) {
    camera = std::make_unique<Camera>(json.at("camera"));
    scene  = std::make_unique<Scene>(json);

    sampler = SamplerFactory::loadSampler(getOptional(json, "sampler", Json()), scene->options.spp, camera->image->resoulation());
    // sampler = std::make_shared < UniformSampler >();

    //load integrator
    {
        const Json&        integratorJson = json.at("integrator");
        const std::string& type           = getOptional(integratorJson, "type", std::string("path_tracer"));
        if (type == "path_tracer") {
            std::string lightSampleStrategy = getOptional(integratorJson, "light_sample", std::string("uniform"));
            integrator                      = std::make_unique<PathIntegrator>(camera, sampler, integratorJson);
        } else if (type == "vol_path_tracer") {
            std::string lightSampleStrategy = getOptional(integratorJson, "light_sample", std::string("uniform"));
            integrator                      = std::make_unique<VolPathIntegrator>(camera, sampler, integratorJson);

        } else if (type == "sppm") {
            int maxBounces = getOptional(integratorJson, "max_bounces", 8);
            integrator     = std::make_unique<PhotonMapper>(camera, integratorJson);
        } else if (type == "bidirectional_path_tracer") {
            integrator = std::make_unique<BDPTIntegrator>(camera, sampler, integratorJson);
        }
    }
}

void Render::Go() {
    TimeCounter counter("Rendering");
    scene->build();
    integrator->process(*scene, *sampler);
    integrator->render(*scene);
    counter.printTimeCount();
}

void Render::renderScene(const std::string path) {
    std::filesystem::path p(path);
    FileUtils::WorkingDir =     p.parent_path().string() + "\\";
    std::ifstream scene_file(path);
    Json json;
    scene_file >> json;
    scene_file.close();
    Render render(json);
    render.Go();
}