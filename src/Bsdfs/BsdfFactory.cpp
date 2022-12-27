#include "Forward.hpp"
#include "BsdfFactory.hpp"
#include "Reflection.hpp"
#include "Dielectric.hpp"
#include "Plastic.hpp"
#include "Hair.h"

#include "Common/Json.hpp"
#include "Texture/TextureFactory.hpp"

#include <sstream>
#include <spdlog/spdlog.h>

namespace BSDFFactory {

    typedef BSDF Material;

    static Spectrum DefaultALbedo = Spectrum(0.9, 0.9, 0.9);


    std::tuple < std::shared_ptr < Texture < Float>>, std::shared_ptr < Texture < Float>>, std::shared_ptr < Texture < Float>> >
    loadRoughness(const Json & json) {
        std::shared_ptr < Texture < Float>> roughness = nullptr, uroughness = nullptr, vroughness = nullptr;
        roughness = TextureFactory::LoadTexture < Float >(json, "roughness", 0.01f);
        uroughness = TextureFactory::LoadTexture < Float >(json, "urounghness");
        vroughness = TextureFactory::LoadTexture < Float >(json, "vrounghness");
        return std::make_tuple(roughness, uroughness, vroughness);
    }


    std::shared_ptr < Material > LoadLambertainMaterial(const Json & j) {
        return std::make_shared < LambertainR >();
    }

    std::shared_ptr < Material > LoadMirrorMaterial(const Json & j) {
        return std::make_shared < SpecularR >();
    }

    std::shared_ptr < Material > LoadDielectricMaterial(const Json & j) {
        bool enalbeT = getOptional(j, "enable_refraction", true);
        Float ior = getOptional(j,"ior",1.33);
        return std::make_shared < Dielectric >(ior, enalbeT);
    }

    std::shared_ptr < Material > LoadConductorMaterial(const Json & j) {
        if ( contains(j, "material") ) {
            return std::make_shared < Conductor >(j["material"]);
        }
        vec3 eta = getOptional(j, "eta", vec3(0.2004376970f, 0.9240334304f, 1.1022119527f));
        vec3 k = getOptional(j, "k", vec3(3.9129485033f, 2.4528477015f, 2.1421879552f));
        return std::make_shared < Conductor >(eta, k);
    }
    std::shared_ptr < Material > LoadPlasticMaterial(const Json & j) {
        std::shared_ptr < Texture < Spectrum>> specularReflection = TextureFactory::LoadTexture < Spectrum >(
                j, "specular_reflection", Spectrum(1));
        std::shared_ptr < Texture < Spectrum>> diffuseReflection = TextureFactory::LoadTexture < Spectrum >(
                j, "diffuse_reflection", Spectrum(0.5));
        Float ior = getOptional(j, "ior", 1.5);
        return std::make_shared <Plastic>(diffuseReflection,specularReflection,ior);
    }

    std::shared_ptr < Material > LoadRoughConductorMaterial(const Json & j) {
        vec3 eta = getOptional(j, "eta", vec3(0.2004376970f, 0.9240334304f, 1.1022119527f));
        vec3 k = getOptional(j, "k", vec3(3.9129485033f, 2.4528477015f, 2.1421879552f));
        std::string distribStr = getOptional(j, "distribution", std::string("beckmann"));
        auto roughnessTuple = loadRoughness(j);
        auto roughCounductorMaterial = std::make_shared < RoughConductor >(eta, k,
                                                                           LoadMicrofacetDistribution(distribStr),
                                                                           std::get < 0 >(roughnessTuple),
                                                                           std::get < 1 >(roughnessTuple),
                                                                           std::get < 2 >(roughnessTuple));
        if ( contains(j, "material") )
            roughCounductorMaterial->setCoundctorByName(j["material"]);

        return roughCounductorMaterial;
    }


//todo
    std::shared_ptr < Material > LoadRoughDielectricMaterial(const Json & j) {
        Float ior = getOptional(j, "ior", 1.5);
        std::string distribStr = getOptional(j, "distribution", std::string("beckmann"));
        auto roughnessTuple = loadRoughness(j);
        auto roughDielectricMaterial = std::make_shared < RoughDielectric >(ior, LoadMicrofacetDistribution(distribStr),
                                                                            std::get < 0 >(roughnessTuple),
                                                                            std::get < 1 >(roughnessTuple),
                                                                            std::get < 2 >(roughnessTuple));

        return roughDielectricMaterial;
    }

    std::shared_ptr < Material > LoadRoughPlasticMaterial(const Json & j) {
        std::shared_ptr < Texture < Spectrum>> specularReflection = TextureFactory::LoadTexture < Spectrum >(
                j, "specular_reflection", Spectrum(1));
        std::shared_ptr < Texture < Spectrum>> diffuseReflection = TextureFactory::LoadTexture < Spectrum >(
                j, "diffuse_reflection", Spectrum(0.5));
        Float ior = getOptional(j, "ior", 1.5);
        std::string distribStr = getOptional(j, "distribution", std::string("beckmann"));
        auto roughnessTuple = loadRoughness(j);
        auto roughPlasticMaterial = std::make_shared < RoughPlastic >(diffuseReflection,specularReflection ,ior,
                                                                      LoadMicrofacetDistribution(distribStr),
                                                                      std::get < 0 >(roughnessTuple),
                                                                      std::get < 1 >(roughnessTuple),
                                                                      std::get < 2 >(roughnessTuple));
        return roughPlasticMaterial;
    }

    std::shared_ptr < Material > LoadHairMaterial(const Json & j ) {
        return std::make_shared <Hair>(j);
    }

    std::shared_ptr < Material > LoadForwardMaterial(const Json & j ) {
        return std::make_shared <ForwardBSDF>();
    }


        std::shared_ptr < Material > LoadDefualtMaterial( ) {
        std::shared_ptr < Material > material = std::make_shared < LambertainR >();
        material->setAlbedo(std::make_shared < ConstantTexture < Spectrum>>(DefaultALbedo));
        return material;
        //   material->Add(new LambertainT(DefaultALbedo));
    }

    std::unordered_map < std::string,
            std::function < std::shared_ptr < Material >(const Json & j)>>
            MaterialLoadTable = {
            {"lambert",          LoadLambertainMaterial},
            {"mirror",           LoadMirrorMaterial},
            {"dielectric",       LoadDielectricMaterial},
            {"conductor",        LoadConductorMaterial},
            {"plastic",        LoadPlasticMaterial},
            {"rough_conductor",  LoadRoughConductorMaterial},
            {"rough_dielectric", LoadRoughDielectricMaterial},
            {"rough_plastic", LoadRoughPlasticMaterial},
            {"hair", LoadHairMaterial},
            {"forward", LoadForwardMaterial}
    };


/**********************************************************************************************************************/
/*********************************************************************************************************************/

    std::shared_ptr < BSDF > LoadBsdfFromJson(const Json & j) {
        // Spectrum  albedo= getOptional(j,"albedo",DefaultALbedo);
        std::shared_ptr < Texture < Spectrum>> albedoTexture = TextureFactory::LoadTexture(j,"albedo",DefaultALbedo);
        std::shared_ptr < Material > material;
        if ( MaterialLoadTable.contains(j["type"]) ) {
            material = MaterialLoadTable[j["type"]](j);
        }
            //todo  support other bsdfs
        else {
            spdlog::info("{} bsdf not loaded correctly.Used Default Bsdf",
                         j["type"]);
            material = LoadDefualtMaterial();
        }

        material->setAlbedo(albedoTexture);

        return material;
    }



    std::unordered_map < std::string, std::shared_ptr < BSDF>>
    LoadBsdfsFromJson(const Json & j) {
        //spdlog::info(to_string(j));
        std::unordered_map < std::string, std::shared_ptr < BSDF>> bsdfMaps;

        for ( auto & bsdfJson: j ) {
            std::string bsdf_name = bsdfJson["name"];
            spdlog::info(bsdf_name);
            auto bsdf = LoadBsdfFromJson(bsdfJson);
            bsdf->name = bsdf_name;
            bsdfMaps[bsdf_name] = bsdf;
        }
        bsdfMaps["default"] = LoadDefualtMaterial();
        spdlog::info(bsdfMaps.size());
        return bsdfMaps;
    }


    std::unordered_map < std::string, std::shared_ptr < BSSRDF>> LoadBssrdfsFromJson(const Json & j) {
        std::unordered_map < std::string, std::shared_ptr < BSSRDF>>  bssrdfMaps;
        for ( auto & bssrdfJson: j ) {
            {
                std::string name = bssrdfJson["name"];
                std::string type  = getOptional(bssrdfJson,"type",std::string("table"));
                Float eta = getOptional(bssrdfJson,"eta",1.5);
                if(type == "disney")
                {
                    auto scatterDistance = TextureFactory::LoadTexture(bssrdfJson,"scatter_distance",Spectrum(0.5));
                    auto color = TextureFactory::LoadTexture(bssrdfJson,"color",Spectrum(0.5));
                    bssrdfMaps[name] = std::make_shared<DisneyBSSRDF>(scatterDistance, color, eta);
                }
                else if(type == "table"){
                    Spectrum sigmas(0.5),sigmaa(0.5);
                    if(bssrdfJson.contains("material"))
                    {
                        GetMediumScatteringProperties(bssrdfJson["material"],&sigmaa,&sigmas);
                    }
                    auto sigmaS = TextureFactory::LoadTexture(bssrdfJson,"sigma_s",sigmas);
                    auto sigmaA = TextureFactory::LoadTexture(bssrdfJson,"sigma_a",sigmaa);

                    Float scale = getOptional(bssrdfJson,"scale",1.0);
                    sigmaS->setScale(scale);sigmaA->setScale(scale);
                    Float g = getOptional(bssrdfJson,"g",0);
                    bssrdfMaps[name] = std::make_shared<TabelBSSRDF>(eta,sigmaA,sigmaS,g);
                }
            }
        }
        return bssrdfMaps;
    }
}