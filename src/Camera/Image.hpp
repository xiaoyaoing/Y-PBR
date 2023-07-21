#pragma  once

#include "Common/Json.hpp"
#include "../Common/math.hpp"
#include "ToneMap.hpp"
#include "spdlog/spdlog.h"

inline vec3 gammaCompress(const vec3 &in) {

    vec3 out;
    for (uint8_t c = 0; c < 3; c++) {
        out[c] = in[c] <= 0.0031308 ? 12.92 * in[c] : 1.055 * std::pow(in[c], 1.0 / 2.4) - 0.055;
    }

//    spdlog::info("{0} {1}", toColorStr(in), toColorStr(out));
    return out;
}

class Image {
    typedef std::atomic<Float> vec3a[3];

    struct Pixel {
        Pixel() {
            Pixel(Spectrum(0));
        }


        Pixel(vec3 s) {
            rgb[0] = s.x;
            rgb[1] = s.y;
            rgb[2] = s.z;
        }

        vec3 add(vec3 value) {
            atomicAdd(rgb[0], value.x);
            atomicAdd(rgb[1], value.y);
            atomicAdd(rgb[2], value.z);
            return vec3();
        }

        vec3 value() const {
            return vec3(rgb[0].load(), rgb[1].load(), rgb[2].load());
        }

        void operator=(vec3 s) {
            rgb[0] = s.x;
            rgb[1] = s.y;
            rgb[2] = s.z;
        }

    protected:
        vec3a rgb;

    };

    vec3 getPixel(int x, int y) const;

    uint32 getIndex(uint32 x, uint32 y) const;

    std::vector<Pixel> buffers;
    std::vector<uint32_t> sampleCounts;
    ToneMap::ToneMapType _tonemapType;

    uint32 _width;
    uint32 _height;
public:
    Image(const ivec2 &res, ToneMap::ToneMapType toneMapType = ToneMap::Filmic)
            : _width(res.x), _height(res.y), _tonemapType(toneMapType), buffers(res.x * res.y),
              sampleCounts(res.x * res.y) {
    }

    Image(const Image &another) {
        _width = another._width;
        _height = another._height;
        _tonemapType = another._tonemapType;
        buffers = std::vector<Pixel>(_width * _height);
        sampleCounts.assign(buffers.size(), 0);
    }

    static void atomicAdd(std::atomic<float> &dst, float add) {
        float current = dst.load();
        float desired = current + add;
        while (!dst.compare_exchange_weak(current, desired))
            desired = current + add;
    }

    void addPixel(uint32 x, uint32 y, vec3 rgb, bool count = false);

    void saveTXT(const std::string &fileName) const;

    void save(const std::string &fileName,Float scale,bool overwrite = false) const;



    ivec2 resoulation() const { return ivec2(_width, _height); }

    int width() const { return _width; }

    int height() const { return _height; }

    inline int product() const { return width() * height(); }

    inline void fill(const Spectrum &spectrum) {
        for (auto &pixel: buffers)
            pixel = vec3();
    }

    vec3 getPixel(int idx) const;


};