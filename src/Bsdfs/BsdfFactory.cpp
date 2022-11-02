#include "BsdfFactory.hpp"
#include "Reflection.hpp"
#include <spdlog/spdlog.h>
#include "../Common/Json.hpp"
#include "Texture/TextureFactory.hpp"
#include "sstream"

namespace BSDFFactory {

typedef   BSDF Material;

static Spectrum  DefaultALbedo =  Spectrum(0.9,0.9,0.9);



std::shared_ptr<Material> LoadLambertainMaterial(Json & j){

    return  std::make_shared<LambertainR>();
   //  material->Add(new LambertainT(albedo));
}

std::shared_ptr<Material> LoadMirrorMaterial(Json j){
    return std::make_shared <SpecularR>();
    //todo support fresnel specular and uRoughness
}

std::shared_ptr<Material> LoadDielectricMaterial(Json j){
    bool enalbeT = getOptional(j,"enable_refraction",true);
    Float ior   = j["ior"];
    return std::make_shared <Dielectric>(ior,enalbeT);
}

std::shared_ptr<Material> LoadConductorMaterial(Json j){
        if(contains(j,"material"))  {
            return std::make_shared <Conductor>(j["material"]);
        }
        vec3 eta  = getOptional(j,"eta",vec3(0.2004376970f, 0.9240334304f, 1.1022119527f));
        vec3 k  =   getOptional(j,"k",vec3(3.9129485033f, 2.4528477015f,2.1421879552f));
        return std::make_shared <Conductor> (eta,k);
}

std::shared_ptr<Material> LoadRoughConductorMaterial(Json j){
    vec3 eta  = getOptional(j,"eta",vec3(0.2004376970f, 0.9240334304f, 1.1022119527f));
    vec3 k  =   getOptional(j,"k",vec3(3.9129485033f, 2.4528477015f,2.1421879552f));
    std::string distribStr = getOptional(j,"distribution",std::string("beckMann"));
    MircroDistributionEnum distribEnum = getOptional("j","distribution",Beckmann);
    std::shared_ptr<Texture<Float>>  roughness = TextureFactory::LoadTexture<Float>(j["roughness"],0.f);
    std::shared_ptr<Texture<Float>>  urounghness =  TextureFactory::LoadTexture<Float>(j["urounghness"]);
    std::shared_ptr<Texture<Float>>  vrounghness =  TextureFactory::LoadTexture<Float>(j["vrounghness"]);

    auto roughCounductorMaterial = std::make_shared<RoughConductor>(eta,k,distribEnum,roughness,urounghness,vrounghness);

    if( contains(j,"material"))
        roughCounductorMaterial->setCoundctorByName(j["material"]);

    return roughCounductorMaterial;
}


//todo
std::shared_ptr<Material> LoadRoughDielectricMaterial(Json j) {
    return nullptr;
}






std::shared_ptr<Material> LoadDefualtMaterial(){
   std::shared_ptr<Material> material =  std::make_shared <LambertainR>();
   material->setAlbedo(std::make_shared<ConstantTexture<Spectrum>>(DefaultALbedo));
   return material;
 //   material->Add(new LambertainT(DefaultALbedo));
}



std::unordered_map<std::string,
            std::function<std::shared_ptr<Material>(Json & j)>>
            MaterialLoadTable  = {
            {"lambert" , LoadLambertainMaterial},
            {"mirror", LoadMirrorMaterial},
            {"dielectric", LoadDielectricMaterial},
            {"conductor", LoadConductorMaterial},
            {"rough_conductor", LoadRoughConductorMaterial},
            {"rough_dielectric", LoadRoughDielectricMaterial}
    };


/**********************************************************************************************************************/
/*********************************************************************************************************************/

std::shared_ptr <BSDF> LoadBsdfFromJson(Json j) {
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
    LoadBsdfsFromJson(Json j) {
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