#include <fstream>
#include "Image.hpp"
#include "iostream"
#include <spdlog/spdlog.h>
#include "../Common/histogram.hpp"
#include "../Common/util.hpp"
#include "IO/ImageIO.hpp"
void Image::savePPM(const std::string &  outPutPath) const {
        std::ofstream file(outPutPath+".ppm");

        file << "P3" << std::endl;
        file << width << " " << height << std::endl;
        file << "255" << std::endl;

        for (unsigned int i = 0; i < height; ++i) {
            for (unsigned int j = 0; j < width; ++j) {
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
        spdlog::info("Write to {0}",outPutPath+".ppm");
}

void Image::saveTXT(const std::string & outPutPath) const {
    std::ofstream file(outPutPath + "txt");
    for ( unsigned int i = 0 ; i < height ; ++ i ) {
        for ( unsigned int j = 0 ; j < width ; ++ j ) {
            const vec3 rgb = getPixel(j, i);
            file << toColorStr(rgb) << std::endl;
        }
        file.close();
    }
}

void Image::saveBMP(const std::string & outPutPath) const {

}


void Image::saveTGA( const std::string &  outPutPath) const {
    HeaderTGA header((uint16_t)width, (uint16_t)height);
    std::ofstream file(outPutPath + ".tga", std::ios::binary);
    file.write(reinterpret_cast<char*>(&header), sizeof(header));
    vec3 averageRadiance(0);
    for (unsigned int i = 0; i < height; ++i) {
        for (unsigned int j = 0; j < width; ++j){
            const vec3 rgb = getPixel(j, i);
            averageRadiance+=rgb/Float(pixels.size());
            vec3 c = glm::clamp(rgb, vec3(0.0), vec3(255.0));
            std::vector<uint8_t> fp={ (uint8_t)c.b, (uint8_t)c.g, (uint8_t)c.r };
            file.write(reinterpret_cast<char*>(fp.data()), fp.size() * sizeof(uint8_t));
        }
    }

    file.close();
    spdlog::info("Average Radiance {0} {1} {2}",averageRadiance.r,averageRadiance.g,averageRadiance.b);
    spdlog::info("Write to {0}",outPutPath+".tga");
}


void Image::savePNG(const std::string &  outPutPath ) const {
    std::vector<unsigned  char>  image;
    image.resize(4*width*height);
    vec3 averageRadiance(0);
    for(unsigned y = 0; y < height; y++)
        for(unsigned x = 0; x < width; x++)
        {
            vec3  rgb = getPixel(x,y);
            averageRadiance+=rgb/Float(width * height);
            if(rgb.x>0)
            {
                int k=1;
            }
            image[4 * width * y + 4 * x + 0] = uint8(rgb.r);
            image[4 * width * y + 4 * x + 1] = uint8(rgb.g);
            image[4 * width * y + 4 * x + 2] = uint8(rgb.b);
            image[4 * width * y + 4 * x + 3] = 255;
        }
    spdlog::info("Average Radiance {0} {1} {2}",averageRadiance.r,averageRadiance.g,averageRadiance.b);
    ImageIO::savePng(outPutPath+".png",image,width,height,4);
    spdlog::info("Write to {0}",outPutPath+".png");
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
    if(rgb.x>0 && abs(rgb.x-0.222587913)>0.000001){
        int k =1;
    }
    if(rgb.x<0 || rgb.y<0 || rgb.z<0 || rgb.x>1|| rgb.y>1 || rgb.z>1){
      //  spdlog::error("Invalid radiance R:{} G:{} B:{}",rgb.x,rgb.y,rgb.z);
//        rgb.x=std::min(rgb.x,Float(1));
//        rgb.y=std::min(rgb.y,Float(1));
//        rgb.z=std::min(rgb.z,Float(1));
    }
    pixels[idx].rgb+=rgb;
}

Image::Image(const ivec2 &  res) {
    width=res.x;
    height=res.y;
    plain=false;
//    Float exposure_EV = getOptional(j, "exposure_compensation", 0.0);
//    Float gain_EV = getOptional(j, "gain_compensation", 0.0);
//    exposure_scale = std::pow(2, exposure_EV);
//    gain_scale = std::pow(2, gain_EV);
    toneMapType=ToneMap::Filmic;
    pixels.resize(width*height,Pixel());
}

void Image::dividePixel(uint32 x, uint32 y, uint32 count) {
    auto t=pixels[getIndex(x,y)].rgb;
    if(t.x>0){
        int k=1;
    }
    pixels[getIndex(x,y)].rgb/=Float(count);
}

void Image::postProgress(){
    Float exposure_factor = plain ? 1.0 : getExposure() * exposure_scale;
    Float gain_factor = plain ? 1.0 : getGain(exposure_factor) * gain_scale;

    spdlog::info("Tone Mapping Type: {0} plain {1}",toneMapType,plain);
    spdlog::info("Exposure {0} Gain {1} ",exposure_factor,gain_factor);
    for(int i=0;i<pixels.size();i++){
        pixels[i].rgb = clamp(255.f * ToneMap::toneMap(toneMapType,pixels[i].rgb),vec3(0),255.f*vec3(1));
      // pixels[i].rgb = 255.f * pixels[i].rgb;
      // pixels[i].rgb  = gammaCompress(ToneMap::toneMap(toneMapType,pixels[i].rgb * exposure_factor) * gain_factor);

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

