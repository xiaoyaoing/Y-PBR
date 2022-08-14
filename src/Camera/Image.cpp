#include <fstream>
#include "Image.hpp"
#include "iostream"
#include <spdlog/spdlog.h>
#include "../Common/histogram.hpp"
#include "../Common/util.hpp"

void Image::savePPM() const {
        std::ofstream file(outputFileName+".ppm");

        file << "P3" << std::endl;
        file << width << " " << height << std::endl;
        file << "255" << std::endl;

        for (unsigned int i = 0; i < height; ++i) {
            for (unsigned int j = 0; j < width; ++j) {
                const vec3 rgb = getPixel(j, i);
                const unsigned int R =
                        std::clamp(static_cast<unsigned int>(255.0f * rgb[0]), 0u, 255u);
                const unsigned int G =
                        std::clamp(static_cast<unsigned int>(255.0f * rgb[1]), 0u, 255u);
                const unsigned int B =
                        std::clamp(static_cast<unsigned int>(255.0f * rgb[2]), 0u, 255u);
                file << R << " " << G << " " << B << std::endl;
            }
        }
        file.close();
        spdlog::info("Write to {0}",outputFileName+".ppm");

}

void Image::saveTGA( ) const {
    HeaderTGA header((uint16_t)width, (uint16_t)height);
    std::ofstream file(outputFileName + ".tga", std::ios::binary);
    file.write(reinterpret_cast<char*>(&header), sizeof(header));
    vec3 averageRadiance;
    for (unsigned int i = 0; i < height; ++i) {
        for (unsigned int j = 0; j < width; ++j){
            const vec3 rgb = getPixel(j, i);
            averageRadiance+=rgb/Float(pixels.size());
            vec3 c = glm::clamp(rgb, vec3(0.0), vec3(1.0)) * Float(255.0);
            std::vector<uint8_t> fp={ (uint8_t)c.b, (uint8_t)c.g, (uint8_t)c.r };
            file.write(reinterpret_cast<char*>(fp.data()), fp.size() * sizeof(uint8_t));
        }
    }

    file.close();
    spdlog::info("Average Radiance {0} {1} {2}",averageRadiance.r,averageRadiance.g,averageRadiance.b);
    spdlog::info("Write to {0}",outputFileName+".tga");

}




vec3 Image::getPixel(int x, int y) const {
    const unsigned int idx = getIndex(x,y);
    return pixels[idx].rgb;
}

uint32 Image::getIndex(uint32 x, uint32 y) const {
   // std::cout<<3 * x + 3 * width * y<<" ";
    return  x +  width * y;
}

void Image::addPixel(uint32 x, uint32 y, vec3 rgb) {
    auto idx = getIndex(x,y);

    if(rgb.x<0 || rgb.y<0 || rgb.z<0 || rgb.x>1|| rgb.y>1 || rgb.z>1){
        spdlog::error("Invalid radiance R:{} G:{} B:{}",rgb.x,rgb.y,rgb.z);
//        rgb.x=std::min(rgb.x,Float(1));
//        rgb.y=std::min(rgb.y,Float(1));
//        rgb.z=std::min(rgb.z,Float(1));
    }

    pixels[idx].rgb=rgb*1.5f;
}

Image::Image(nlohmann::json j) {
    width=j.at("width");
    height=j.at("height");
    outputFileName=j.at("outputFileName");
    plain=false;

    Float exposure_EV = getOptional(j, "exposure_compensation", 0.0);
    Float gain_EV = getOptional(j, "gain_compensation", 0.0);
    exposure_scale = std::pow(2, exposure_EV);
    gain_scale = std::pow(2, gain_EV);

    toneMapType=ToneMap::Aces;
    pixels.resize(width*height);
}

void Image::dividePixel(uint32 x, uint32 y, uint32 count) {
    pixels[getIndex(x,y)].rgb/=Float(count);
}

void Image::postProgress(){
    Float exposure_factor = plain ? 1.0 : getExposure() * exposure_scale;
    Float gain_factor = plain ? 1.0 : getGain(exposure_factor) * gain_scale;

    spdlog::info("Tone Mapping Type: {0}",toneMapType);
    spdlog::info("Exposure {0} Gain {1}",exposure_factor,gain_factor);
    for(int i=0;i<pixels.size();i++){
        pixels[i].rgb=gammaCompress(
                            ToneMap::toneMap(toneMapType,pixels[i].rgb * exposure_factor) * gain_factor);
       // pixels[i].rgb=ToneMap::toneMap(ToneMap::Filmic,pixels[i].rgb);
    }
}

/*******************************************************************************************
Histogram method to find the intensity level L that 50% of the pixels has higher/lower intensity than.
The returned exposure factor is then 0.5/L, which if multiplied by each pixel in the image will make
50% of the pixels be < 0.5 and 50% of the pixels be > 0.5.
********************************************************************************************/
Float Image::getExposure() const
{
    std::vector<Float> brightness(pixels.size());
    for (size_t i = 0; i < pixels.size(); i++)
    {
        brightness[i] = compAdd(pixels[i].rgb) / 3.0;
    }
    Histogram histogram(brightness, 65536);
    Float L = histogram.level(0.5);
    return L > 0.0 ? 0.5 / L : 1.0;
}

/**************************************************************************
Histogram method to find the gain that positions the histogram to the right
***************************************************************************/
Float Image::getGain(Float exposure_factor) const
{
    std::vector<Float> brightness(pixels.size());
    for (size_t i = 0; i < pixels.size(); i++)
    {
        brightness[i] = compAdd(ToneMap::toneMap
                                        (toneMapType,pixels[i].rgb * exposure_factor)) / Float(3.0);
    }
    Histogram histogram(brightness, 65536);
    Float L = histogram.level(0.99);
    return L > 0.0 ? 0.99 / L : 1.0;
}

