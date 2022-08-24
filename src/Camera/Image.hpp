#pragma  once

#include <nlohmann/json.hpp>
#include "../Common/math.hpp"
#include "ToneMap.hpp"
#include "spdlog/spdlog.h"

inline  vec3 gammaCompress(const vec3 &in)
{

    vec3 out;
    for (uint8_t c = 0; c < 3; c++)
    {
        out[c] = in[c] <= 0.0031308 ? 12.92 * in[c] : 1.055 * std::pow(in[c], 1.0 / 2.4) - 0.055;
    }

//    spdlog::info("{0} {1}", toColorStr(in), toColorStr(out));
    return out;
}

class Image{

    struct Pixel{
        vec3 rgb = vec3(0);
    };

    /**************************************************************************
    Hard coded (except for dimensions) uncompressed 24bpp true-color TGA header.
    After writing this to file, the RGB bytes can be dumped in sequence
    (left to right, top to bottom) to create a TGA image.
    ***************************************************************************/
    struct HeaderTGA
    {
        HeaderTGA(uint16_t width, uint16_t height)
                : width(width), height(height) {}

    private:
        uint8_t begin[12] = { 0, 0, 2 };
        uint16_t width;
        uint16_t height;
        uint8_t end[2] = { 24, 32 };
    };

    vec3 getPixel(int x,int y) const ;

    uint32 getIndex(uint32 x,uint32 y) const ;

    std::string  outputFileName;
    std::vector<Pixel> pixels;
    ToneMap::ToneMapType toneMapType;
    bool plain;
    Float exposure_scale;
    Float gain_scale;


public:
    Image(nlohmann::json json);

    void addPixel(uint32 x,uint32 y,vec3 rgb)  ;

    void dividePixel(uint32 x,uint32 y,uint32);

    void savePPM() const;

    void saveTGA() const;

    void savePNG() const ;

    void postProgress();

    uint32 width;
    uint32 height;


    Float getExposure( ) const;

    Float getGain(Float exposure_factor) const;
};