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
    //    for (int i = 0; i < 20; i++) {
    //        for (int j = 0; j < 20; j++) {
    //            Float beta_m = 0.01 + i * 0.01 * (i / 4 +1) * i;
    //            Float beta_n = 0.01 + j * 0.01 * (j / 4 + 1) * j;
    //            beta_m = beta_n;
    //
    //            std::cout << "beta_m";
    FileUtils::WorkingDir = path;
    std::ifstream scene_file(FileUtils::WorkingDir + "scene.json");

    nlohmann::json json;
    scene_file >> json;
    scene_file.close();
    //  auto t = json["bsdfs"][0];

    //            json["bsdfs"][0]["beta_m"] = beta_m;
    //            json["bsdfs"][0]["beta_n"] = beta_n;
    //            char s[100];
    //            sprintf_s(s, sizeof(s), "betaM%.2f betaN%.2fTT.png", beta_m, beta_n);

    //   json["renderer"]["output_file"] = std::string(s);

    Render render(json);

    render.Go();
    // }
}