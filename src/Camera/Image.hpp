#pragma once

#include "Common/Json.hpp"
#include "../Common/math.hpp"
#include "ToneMap.hpp"
#include "spdlog/spdlog.h"

inline vec3 gammaCompress(const vec3& in) {

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

        Pixel(const Pixel& another) {
            rgb[0] = another.rgb[0].load();
            rgb[1] = another.rgb[1].load();
            rgb[2] = another.rgb[2].load();
        }

        Pixel& operator=(const Pixel& another) {
            if (this == &another) {
                return *this;// Handle self-assignment
            }
            rgb[0] = another.rgb[0].load();
            rgb[1] = another.rgb[1].load();
            rgb[2] = another.rgb[2].load();
            return *this;
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
            if (isnan(rgb[0].load())) {
                spdlog::info("error");
            }
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

    uint32 getIndex(uint32 x, uint32 y) const;

    std::vector<Pixel>    buffers;
    std::vector<uint32_t> sampleCounts;
    ToneMap::ToneMapType  _tonemapType;

    uint32 _width;
    uint32 _height;

public:
    Image(const ivec2& res, ToneMap::ToneMapType toneMapType = ToneMap::Filmic)
        : _width(res.x), _height(res.y), _tonemapType(toneMapType), buffers(res.x * res.y, Pixel()),
          sampleCounts(res.x * res.y, 0) {
    }

    Image(const Image& another) {
        _width       = another._width;
        _height      = another._height;
        _tonemapType = another._tonemapType;
        buffers      = std::vector<Pixel>(_width * _height);
        sampleCounts.assign(buffers.size(), 0);
    }

    static void atomicAdd(std::atomic<float>& dst, float add) {
        float current = dst.load();
        float desired = current + add;
        while (!dst.compare_exchange_weak(current, desired))
            desired = current + add;
    }

    void addPixel(uint32 x, uint32 y, vec3 rgb, bool count = true);

    void saveTXT(const std::string& fileName) const;

    void save(const std::string& fileName, Float scale, bool overwrite = false) const;

    void linerarNormalize() {
        // return ;
        vec3 minVal(1e5f), maxVal(-1e5f);
        for (int x = 0; x < width(); x++) {
            for (int y = 0; y < height(); y++) {
                int idx = getIndex(x, y);
                if (sampleCounts[idx]) {
                    minVal = min(minVal, buffers[idx].value() / Float(sampleCounts[idx]));
                    maxVal = max(maxVal, buffers[idx].value() / Float(sampleCounts[idx]));
                }
            }
        }
        auto diffVal = maxVal - minVal;
        for (int i = 0; i < buffers.size(); i++) {
            if (sampleCounts[i] <= 0)
                continue;
            auto value = getPixel(i);
            auto temp  = value;
            value      = (value - minVal) / diffVal;
            buffers[i] = Pixel(value * Float(sampleCounts[i]));
        }
    }

    void normalize() {
        return;
        vec3 average(0);
        int  count = 0;
        for (int x = 0; x < width(); x++) {
            for (int y = 0; y < height(); y++) {
                int idx = getIndex(x, y);
                if (sampleCounts[idx]) {
                    average = average * Float(count) / Float(count + 1) + buffers[idx].value() / Float(count + 1);
                    count++;
                }
            }
        }
        vec3 diff_squared_i;
        count = 0;
        for (int x = 0; x < width(); x++) {
            for (int y = 0; y < height(); y++) {
                int idx = getIndex(x, y);
                if (sampleCounts[idx]) {
                    diff_squared_i *= Float(count) / Float(count + 1);
                    diff_squared_i += sqr(buffers[idx].value() - average) / Float(count + 1);
                    count++;
                }
            }
        }
        diff_squared_i /= count;
        diff_squared_i       = sqrt(diff_squared_i);
        auto normalize_pixel = [&](const Pixel& pixel) {
            auto value = pixel.value();
            value      = (value - average) / diff_squared_i;
            return Pixel(value);
        };
        std::transform(buffers.begin(), buffers.end(), buffers.begin(), normalize_pixel);
    }

    void clear() {
        std::fill(buffers.begin(), buffers.end(), Pixel());
        std::fill(sampleCounts.begin(), sampleCounts.end(), 0);
    }

    ivec2 resoulation() const { return ivec2(_width, _height); }

    int width() const { return _width; }

    int height() const { return _height; }

    inline int product() const { return width() * height(); }

    inline void fill(const Spectrum& spectrum) {
        std::fill(buffers.begin(), buffers.end(), spectrum);
    }

    vec3 getPixel(int idx) const;

    vec3 getPixel(int x, int y) const;
};