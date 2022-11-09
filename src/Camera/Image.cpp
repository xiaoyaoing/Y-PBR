#include "Image.hpp"
#include "Common/histogram.hpp"
#include "Common/Json.hpp"
#include "IO/ImageIO.hpp"
#include "IO/FileUtils.hpp"

#include <spdlog/spdlog.h>
#include <iostream>
#include <fstream>

void Image::savePPM() const {
        std::ofstream file(fileName+".ppm");

        file << "P3" << std::endl;
        file << _width << " " << _height << std::endl;
        file << "255" << std::endl;

        for ( unsigned int i = 0; i < _height; ++i) {
            for ( unsigned int j = 0; j < _width; ++j) {
                const vec3 rgb = getPixel(j, i);
                const unsigned int R =
                        std::clamp(static_cast<unsigned int>( rgb[0]), 0u, 255u);
                const unsigned int G =
                        std::clamp(static_cast<unsigned int>(rgb[1]), 0u, 255u);
                const unsigned int B =
                        std::clamp(static_cast<unsigned int>(rgb[2]), 0u, 255u);
                file << R << " " << G << " " << B << std::endl;
            }
        }
        file.close();
        spdlog::info("Write to {0}",fileName+".ppm");
}

void Image::saveTXT() const {
    std::ofstream file(fileName + "txt");
    for ( unsigned int i = 0 ; i < _height ; ++ i ) {
        for ( unsigned int j = 0 ; j < _width ; ++ j ) {
            const vec3 rgb = getPixel(j, i);
            file << toColorStr(rgb) << std::endl;
        }
        file.close();
    }
}

void Image::saveBMP() const {

}


void Image::saveTGA( ) const {
    HeaderTGA header((uint16_t)_width, (uint16_t)_height);
    std::ofstream file(fileName + ".tga", std::ios::binary);
    file.write(reinterpret_cast<char*>(&header), sizeof(header));
    vec3 averageRadiance(0);
    for ( unsigned int i = 0; i < _height; ++i) {
        for ( unsigned int j = 0; j < _width; ++j){
            const vec3 rgb = getPixel(j, i);
            averageRadiance+=rgb/Float(pixels.size());
            vec3 c = glm::clamp(rgb, vec3(0.0), vec3(255.0));
            std::vector<uint8_t> fp={ (uint8_t)c.b, (uint8_t)c.g, (uint8_t)c.r };
            file.write(reinterpret_cast<char*>(fp.data()), fp.size() * sizeof(uint8_t));
        }
    }

    file.close();
//    spdlog::info("Average Radiance {0} {1} {2}",averageRadiance.r,averageRadiance.g,averageRadiance.b);
//    spdlog::info("Write to {0}",fileName+".tga");
}


void Image::savePNG( ) const {
    std::vector<unsigned  char>  image;
    image.resize(4 * _width * _height);
    vec3 averageRadiance(0);
    for( unsigned y = 0; y < _height; y++)
        for( unsigned x = 0; x < _width; x++)
        {
            vec3  rgb = getPixel(x,y);
            averageRadiance+=rgb/Float(_width * _height);
            image[4 * _width * y + 4 * x + 0] = uint8(rgb.r);
            image[4 * _width * y + 4 * x + 1] = uint8(rgb.g);
            image[4 * _width * y + 4 * x + 2] = uint8(rgb.b);
            image[4 * _width * y + 4 * x + 3] = 255;
        }
    const std::string imagePath = FileUtils::getFilePath(FileUtils::WorkingDir+fileName,"png",false);
    ImageIO::savePng(imagePath, image, _width, _height, 4);
        spdlog::info("Average Radiance {0} {1} {2}",averageRadiance.r,averageRadiance.g,averageRadiance.b);
    //    spdlog::info("Write to {0}",imagePath);
}


vec3 Image::getPixel(int x, int y) const {
    const unsigned int idx = getIndex(x,y);
    return pixels[idx].rgb;
}

uint32 Image::getIndex(uint32 x, uint32 y) const {
   // std::cout<<3 * x + 3 * width * y<<" ";
    return x + _width * y;
}

void Image::addPixel(uint32 x, uint32 y, vec3 rgb) {
    uint32 idx = getIndex(x,y);
    pixels[idx].rgb+=rgb;
}

void Image::dividePixel(uint32 x, uint32 y, uint32 count) {
    pixels[getIndex(x,y)].rgb/=Float(count);
}

void Image::postProgress(){
    for(int i=0;i<pixels.size();i++){
        pixels[i].rgb = clamp(255.f * ToneMap::toneMap(_tonemapType, pixels[i].rgb), vec3(0), 255.f * vec3(1));
    }
}

///*******************************************************************************************
//Histogram method to find the intensity level L that 50% of the pixels has higher/lower intensity than.
//The returned exposure factor is then 0.5/L, which if multiplied by each pixel in the image will make
//50% of the pixels be < 0.5 and 50% of the pixels be > 0.5.
//********************************************************************************************/
//Float Image::getExposure() const
//{
//    std::vector<Float> brightness(pixels.size());
//    for (size_t i = 0; i < pixels.size(); i++)
//    {
//        brightness[i] = compAdd(pixels[i].rgb) / 3.0;
//    }
//
//    Histogram histogram(brightness, 65536);
//    Float L = histogram.level(0.5);
//    return L > 0.0 ? 0.5 / L : 1.0;
//}
//
///**************************************************************************
//Histogram method to find the gain that positions the histogram to the right
//***************************************************************************/
//Float Image::getGain(Float exposure_factor) const
//{
//    std::vector<Float> brightness(pixels.size());
//    for (size_t i = 0; i < pixels.size(); i++)
//    {
//        brightness[i] = compAdd(ToneMap::toneMap
//                                        (toneMapType,pixels[i].rgb * exposure_factor)) / Float(3.0);
//    }
//    Histogram histogram(brightness, 65536);
//    Float L = histogram.level(0.99);
//    return L > 0.0 ? 0.99 / L : 1.0;
//}

