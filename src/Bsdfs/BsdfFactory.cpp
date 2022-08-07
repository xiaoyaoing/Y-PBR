#include "BsdfFactory.hpp"
#include "Lambertain.hpp"
#include "spdlog/spdlog.h"

namespace BsdfFactory {

static Spectrum  DefaultALbedo =  Spectrum(0.5,0.5,0.5);

    std::shared_ptr<Bsdf> getDefaultBsdf(){
    return std::make_shared<LambertainBsdf>(DefaultALbedo);
}

std::shared_ptr < Bsdf > LoadBsdfFromJson(nlohmann::json j) {
    if(j["type"]=="lambert"){
        return CreateLambertainBsdf(j);
    }
    //todo  support other bsdfs
    else{
        spdlog::info("{} bsdf not loaded correctly.Used Default Bsdf",
                     j["type"]);
        return std::make_shared <LambertainBsdf>(DefaultALbedo);

    }

    }

std::unordered_map < std::string, std::shared_ptr< Bsdf>>
    LoadBsdfsFromJson(nlohmann::json j) {
        //spdlog::info(to_string(j));
        std::unordered_map < std::string, std::shared_ptr< Bsdf>> bsdf_maps;
        for(auto bsdf_json:j){
            bsdf_maps[bsdf_json["name"]] = LoadBsdfFromJson(bsdf_json);
        }

        bsdf_maps["default"] =   std::make_shared <LambertainBsdf>(DefaultALbedo);

        return bsdf_maps;
    }
}