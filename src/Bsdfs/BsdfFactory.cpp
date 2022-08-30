#include "BsdfFactory.hpp"
#include "Reflection.hpp"
#include <spdlog/spdlog.h>
#include "../Common/util.hpp"
#include "sstream"
namespace BsdfFactory {


typedef  Bsdf Material;
static Spectrum  DefaultALbedo =  Spectrum(0.9,0.9,0.9);



std::shared_ptr<Material> LoadLambertainMaterial(nlohmann::json & j){
    std::shared_ptr<Material> material =std::make_shared <Material>();
    Spectrum albedo;
    nlohmann::json & r =j["albedo"];
    if (r.type() == nlohmann::json::value_t::string)
    {
        std::string hex_string = r.get<std::string>();
        if (hex_string.size() == 7 && hex_string[0] == '#')
        {
            hex_string.erase(0, 1);
            std::stringstream ss;
            ss << std::hex << hex_string;
            uint32_t color_int;
            ss >> color_int;
            albedo = intToColor(color_int);
        }
    }
    else
    {
       albedo= r.get<Spectrum>();
    }
    material->Add(new LambertainR(albedo));
   //  material->Add(new LambertainT(albedo));
    return material;
}

std::shared_ptr<Material> LoadMirrorMaterial(nlohmann::json j){
    std::shared_ptr<Material> material =std::make_shared <Material>();
    //todo support fresnel specular and uRoughness
    material->Add(new SpecularR());
    return material;
}

std::shared_ptr<Material> LoadDielectricMaterial(nlohmann::json j){
    std::shared_ptr<Material> material =std::make_shared <Material>();
    bool enalbeT = getOptional(j,"enable_refraction",true);
    Float ior   = j["ior"];
    Spectrum  albedo = getOptional(j,"abledo",Spectrum(1));
    material->Add(new Dielectric(ior,albedo,enalbeT));
    return material;
}

std::shared_ptr<Material> LoadConductMaterial(nlohmann::json j){
        std::shared_ptr<Material> material =std::make_shared <Material>();
        Float k  =j["k"];
        Float eta   = j["eta"];
        Spectrum  albedo = getOptional(j,"abledo",Spectrum(1));
        material->Add(new Conductor(albedo,eta,k));
        return material;
    }





std::shared_ptr<Material> LoadDefualtMaterial(){
    std::shared_ptr<Material> material =std::make_shared <Material>();
    material->Add(new LambertainR(DefaultALbedo));
 //   material->Add(new LambertainT(DefaultALbedo));
    return material;
}



std::unordered_map<std::string,
            std::function<std::shared_ptr<Material>(nlohmann::json & j)>>
            MaterialLoadTable  = {
            {"lambert" , LoadLambertainMaterial},
            {"mirror", LoadMirrorMaterial},
            {"dielectric", LoadDielectricMaterial},
            {"conductor", LoadConductMaterial}
    };


/**********************************************************************************************************************/
/*********************************************************************************************************************/

std::shared_ptr <Bsdf> LoadBsdfFromJson(nlohmann::json j) {

    if(MaterialLoadTable.contains(j["type"])){
        return MaterialLoadTable[j["type"]](j);
    }
    //todo  support other bsdfs
    else{
        spdlog::info("{} bsdf not loaded correctly.Used Default Bsdf",
                     j["type"]);
        return LoadDefualtMaterial();
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

        bsdf_maps["default"] =   LoadDefualtMaterial();

        spdlog::info(bsdf_maps.size());
        return bsdf_maps;
    }
}