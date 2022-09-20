#include "BsdfFactory.hpp"
#include "Reflection.hpp"
#include <spdlog/spdlog.h>
#include "../Common/util.hpp"
#include "Texture/TextureFactory.hpp"
#include "sstream"

namespace BSDFFactory {

typedef   BSDF Material;

static Spectrum  DefaultALbedo =  Spectrum(0.9,0.9,0.9);



std::shared_ptr<Material> LoadLambertainMaterial(nlohmann::json & j){

    return  std::make_shared<LambertainR>();
   //  material->Add(new LambertainT(albedo));
}

std::shared_ptr<Material> LoadMirrorMaterial(nlohmann::json j){
    return std::make_shared <SpecularR>();
    //todo support fresnel specular and uRoughness
}

std::shared_ptr<Material> LoadDielectricMaterial(nlohmann::json j){
    bool enalbeT = getOptional(j,"enable_refraction",true);
    Float ior   = j["ior"];
    return std::make_shared <Dielectric>(ior,enalbeT);
}

std::shared_ptr<Material> LoadConductorMaterial(nlohmann::json j){
        if( contains(j,"material"))  {
            return std::make_shared <Conductor>(j["material"]);
        }
        vec3 k  =j["k"];
        vec3 eta   = j["eta"];
        return std::make_shared <Conductor> (eta,k);
}





std::shared_ptr<Material> LoadDefualtMaterial(){
   std::shared_ptr<Material> material =  std::make_shared <LambertainR>();
   material->setAlbedo(std::make_shared<ConstantTexture<Spectrum>>(DefaultALbedo));
   return material;
 //   material->Add(new LambertainT(DefaultALbedo));
}



std::unordered_map<std::string,
            std::function<std::shared_ptr<Material>(nlohmann::json & j)>>
            MaterialLoadTable  = {
            {"lambert" , LoadLambertainMaterial},
            {"mirror", LoadMirrorMaterial},
            {"dielectric", LoadDielectricMaterial},
            {"conductor", LoadConductorMaterial}
    };


/**********************************************************************************************************************/
/*********************************************************************************************************************/

std::shared_ptr <BSDF> LoadBsdfFromJson(nlohmann::json j) {
   // Spectrum  albedo= getOptional(j,"albedo",DefaultALbedo);
    std::shared_ptr<Texture<Spectrum>> albedoTexture = TextureFactory::LoadTexture(j["albedo"],DefaultALbedo);
    //    std::shared_ptr<ConstantTexture<Spectrum>> albedoTexture = std::make_shared<ConstantTexture<Spectrum>>(albedo);
   // std::shared_ptr<Texture<Float>>   bumpMapTexture = std::make_shared<ConstantTexture<Spectrum>>(albedo);
    std::shared_ptr<Material> material ;
    if(MaterialLoadTable.contains(j["type"])){
            material =  MaterialLoadTable[j["type"]](j);
    }
    //todo  support other bsdfs
    else{
        spdlog::info("{} bsdf not loaded correctly.Used Default Bsdf",
                     j["type"]);
        material = LoadDefualtMaterial();
    }

    material->setAlbedo(albedoTexture);

   return material;
}

std::unordered_map < std::string, std::shared_ptr<BSDF>>
    LoadBsdfsFromJson(nlohmann::json j) {
        //spdlog::info(to_string(j));
        std::unordered_map < std::string, std::shared_ptr<BSDF>> bsdf_maps;


        for(auto bsdf_json:j){
            std::string bsdf_name=bsdf_json["name"];
            spdlog::info(bsdf_name);
            auto bsdf=LoadBsdfFromJson(bsdf_json);
            bsdf->name=bsdf_name;
            bsdf_maps[bsdf_name] = bsdf;
        }
        bsdf_maps["default"] =   LoadDefualtMaterial();
        spdlog::info(bsdf_maps.size());
        return bsdf_maps;
    }
}