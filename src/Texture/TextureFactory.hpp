#include "Common/Texture.hpp"
#include "ConstantTexture.hpp"
#include "CheckTexture.hpp"
#include "BitMapTexture.hpp"
#include "Colors/Spectrum.hpp"
#include "Common/math.hpp"
#include "Common/Json.hpp"
#include "Common/Json.hpp"
#include "IO/FileUtils.hpp"
namespace  TextureFactory{
    template <class T>
    std::shared_ptr<Texture<T>> LoadTexture(const Json & texturJson){
        if(texturJson.is_null()){
            return nullptr;
        }
        if(texturJson.is_number()){
            return std::make_shared <ConstantTexture<T>>(texturJson);
        }
        if(texturJson.is_array()){
            return std::make_shared <ConstantTexture<T>>(texturJson);
        }
        if(texturJson.is_string()){
              auto bmt =  std::make_shared <BitMapTexture<T>>(FileUtils::WorkingDir+texturJson.get<std::string>());
              bmt->LoadResources();
              return bmt;
        }
        if(texturJson.is_object()){
            std::string textureType = texturJson.at("type");
            if(textureType == "checker"){
                T onColor= getOptional(texturJson,"on_color",T());
                T offColor= getOptional(texturJson,"off_color",T());
                int resU = getOptional(texturJson,"res_u",20);
                int resV = getOptional(texturJson,"res_v",20);

                return std::make_shared <CheckTexture<T>>(onColor,offColor,resU,resV);
            }
        }
        return nullptr;
       // return std::make_shared <ConstantTexture<T>>(defaultValue);
    }

    template <class T>
    std::shared_ptr<Texture<T>> LoadTexture(const Json & texturJson,T defaultValue){
        if(auto texture = LoadTexture<T>(texturJson))
            return texture;
        return std::make_shared<ConstantTexture<T>>(defaultValue);
    }

    //judge whether the field exists
    template <class T>
    std::shared_ptr<Texture<T>> LoadTexture(const Json & json,const std::string & field){
        if(!json.contains(field))
            return nullptr;
        return LoadTexture <T>(json[field]);
    }

    template <class T>
    std::shared_ptr<Texture<T>> LoadTexture(const Json & json,const std::string & field,T defaultValue){
        if(auto texture = LoadTexture<T>(json,field))
            return  texture;
        return std::make_shared<ConstantTexture<T>>(defaultValue);
    }
}



