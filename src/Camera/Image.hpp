#pragma  once

#include "Common/Json.hpp"
#include "../Common/math.hpp"
#include "ToneMap.hpp"
#include "spdlog/spdlog.h"

inline vec3 gammaCompress(const vec3 & in) {

    vec3 out;
    for ( uint8_t c = 0 ; c < 3 ; c ++ ) {
        out[c] = in[c] <= 0.0031308 ? 12.92 * in[c] : 1.055 * std::pow(in[c], 1.0 / 2.4) - 0.055;
    }

//    spdlog::info("{0} {1}", toColorStr(in), toColorStr(out));
    return out;
}

class Image {

    struct Pixel {
        vec3 rgb = vec3(0);
    };

    /**************************************************************************
    Hard coded (except for dimensions) uncompressed 24bpp true-color TGA header.
    After writing this to file, the RGB bytes can be dumped in sequence
    (left to right, top to bottom) to create a TGA image.
    ***************************************************************************/
    struct HeaderTGA {
        HeaderTGA(uint16_t width, uint16_t height)
                : width(width), height(height) {}

    private:
        uint8_t begin[12] = {0, 0, 2};
        uint16_t width;
        uint16_t height;
        uint8_t end[2] = {24, 32};
    };

    vec3 getPixel(int x, int y) const;

    uint32 getIndex(uint32 x, uint32 y) const;

    std::vector < Pixel > pixels;
    ToneMap::ToneMapType _tonemapType;
    Float exposure_scale;
    Float gain_scale;
    
    std::string  fileName;
    uint32 _width;
    uint32 _height;
public:
    Image(const ivec2 & res,const std::string & fileName,ToneMap::ToneMapType toneMapType = ToneMap::Filmic )
          : _width(res.x), _height(res.y), fileName(fileName), _tonemapType(toneMapType){
        pixels.resize(_width * _height, Pixel());
    }

    void addPixel(uint32 x, uint32 y, vec3 rgb);
    void dividePixel(uint32 x, uint32 y, uint32);


    void savePPM() const;
    void saveTGA() const;
    void savePNG() const;
    void saveBMP() const;
    void saveTXT() const;


    void postProgress( );
    void fill(const Spectrum & spectrum){
        std::fill(pixels.begin(),pixels.end(),Pixel{spectrum});
    }
//    Float getExposure( ) const;
//    Float getGain(Float exposure_factor) const;

    ivec2 resoulation(){ return ivec2(_width, _height);}
    int width(){return _width;}
    int height(){return _height;}
};