#include "BsdfFactory.hpp"
#include "Reflection.hpp"
#include "spdlog/spdlog.h"

namespace BsdfFactory {


typedef  Bsdf Material;


std::shared_ptr<Material> LoadLambertainBsdf(nlohmann::json j){
    std::shared_ptr<Material> material =std::make_shared <Material>();

    Spectrum albedo = j["albedo"];
    material->Add(new Lambertain(albedo));
}


static Spectrum  DefaultALbedo =  Spectrum(0.5,0.5,0.5);

    std::shared_ptr<Bsdf> getDefaultBsdf(){
    return std::make_shared<LambertainBsdf>(DefaultALbedo);
}

std::shared_ptr <Bsdf> LoadBsdfFromJson(nlohmann::json j) {
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
            std::string bsdf_name=bsdf_json["name"];
            auto bsdf=LoadBsdfFromJson(bsdf_json);
            bsdf->name=bsdf_name;
            bsdf_maps[bsdf_name] = bsdf;
        }

        bsdf_maps["default"] =   std::make_shared <LambertainBsdf>(DefaultALbedo);

        spdlog::info(bsdf_maps.size());
        return bsdf_maps;
    }
}