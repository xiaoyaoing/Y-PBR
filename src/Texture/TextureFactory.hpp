#include "Common/Texture.hpp"
#include "ConstantTexture.hpp"
#include "CheckTexture.hpp"
#include "BitMapTexture.hpp"
#include "Colors/Spectrum.hpp"
#include "Common/math.hpp"
#include "nlohmann/json.hpp"
#include "Common/util.hpp"

namespace  TextureFactory{
    template <class T>
    std::shared_ptr<Texture<T>> LoadTexture(const nlohmann::json & texturJson){
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
//            std::string hex_string = texturJson.get<std::string>();
//            if (hex_string.size() == 7 && hex_string[0] == '#')
//            {
//                hex_string.erase(0, 1);
//                std::stringstream ss;
//                ss << std::hex << hex_string;
//                uint32_t color_int;
//                ss >> color_int;
//                Spectrum  albedo = intToColor(color_int);
//                return std::make_shared <ConstantTexture<T>>(albedo);
//            }

              return nullptr;
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
    std::shared_ptr<Texture<T>> LoadTexture(const nlohmann::json & texturJson,T defaultValue){
        if(auto texture = LoadTexture<T>(texturJson))
            return texture;

        return std::make_shared<ConstantTexture<T>>(defaultValue);
    }
}



