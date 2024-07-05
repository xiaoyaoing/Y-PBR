#pragma once

#include "Common/Texture.hpp"
#include "Sampler/Distrib.hpp"
#include "Common/Spectrum.hpp"
#include "IO/ImageIO.hpp"
#include "IO/FileUtils.hpp"
#include  "Common/Log.h"
struct Rgba {
    uint8 c[4];

    vec3 normalize() const {
        return vec3(float(c[0]), float(c[1]), float(c[2])) * (1.0f / 255.0f);
    }
};

template<class T>
class BitMapTexture : public Texture<T> {
public:
    enum class TexelType : uint32 {
        SCALAR_LDR = 0,
        SCALAR_HDR = 1,
        RGB_LDR    = 2,
        RGB_HDR    = 3,
    };

    BitMapTexture(const std::string& path) : _path(path) {
    }

    BitMapTexture(void* texels, int w, int h, TexelType texelType, bool linear, bool clamp) : _texels(texels), _w(w), _h(h), _texelType(texelType), _linear(linear), _clamp(clamp) {
        _isRGB = int(texelType) & 2;
        _isHdr = int(texelType) & 1;
        for (int y = 0; y < _h; ++y)
            for (int x = 0; x < _w; x++) {
                _average += getValue(x, y) / (float)(_w * _h);
            }
    }

    T eval(const Intersection* si) const override {
        return eval(si->uv);
    }

    void makeSamplable(TextureMapJacobian jacobian) override {
        if (_distribution[jacobian])
            return;

        std::vector<float> weights(_w * _h);

        for (int y = 0, idx = 0; y < _h; ++y) {
            float rowWeight = 1.0f;
            if (jacobian == MAP_SPHERICAL)
                rowWeight *= std::sin((y * Constant::PI) / _h);
            for (int x = 0; x < _w; ++x, ++idx)
                weights[idx] = weight(x, y) * rowWeight;
        }
        _distribution[jacobian].reset(new Distribution2D(weights.data(), _w, _h));
    }

    vec2 sample(TextureMapJacobian jacobian, const vec2& uv, Float* pdf) const override {
        vec2 newUv = _distribution[jacobian]->SampleContinuous(uv, pdf);
        *pdf       = *pdf * _w * _h;
        return vec2(newUv.x, 1 - newUv.y);
    }

    Float pdf(TextureMapJacobian jacobian, const vec2& uv) const override {
        vec2  newuv(uv);
        Float res = _distribution[jacobian]->Pdf(vec2(newuv.x, 1 - newuv.y));
        return res;
    }

    void convert() {

        //    ImageIO::saveLdr("a.png", ldr.get(), width(), height(), 3, overwrite);
    }

    T evalWithStochasticSampling(vec2 uv,float u){
        uv[0] = std::max(0.f,uv[0] * width() - 0.5f);
        uv[1] = std::max(0.f,uv[1] * height() - 0.5f);


        int u0 = std::floor(uv[0]), v0 = std::floor(uv[1]);
        Float du = uv[0] - u0, dv = uv[1] - v0;
        if (u < du) {
            ++u0;
            u /= du;
        } else
            u = (u - du) / (1 - du);

        if (u < dv)
            ++v0;
        return getValue(u0, v0);
    }
    T eval(vec2 uv) const override {
        float u = uv.x * _w;
        float v = (1 - uv.y) * _h;
        if (_linear) {
            u -= 0.5f;
            v -= 0.5f;
        }
        int iu0 = u < 0.0f ? -int(-u) - 1 : int(u);
        int iv0 = v < 0.0f ? -int(-v) - 1 : int(v);
        int iu1 = iu0 + 1;
        int iv1 = iv0 + 1;

        u -= iu0;
        v -= iv0;
        if (_clamp) {
            iu0 = clamp(iu0, 0, _w - 1);
            iu1 = clamp(iu1, 0, _w - 1);
            iv0 = clamp(iv0, 0, _h - 1);
            iv1 = clamp(iv1, 0, _h - 1);
        } else {
            iu0 = ((iu0 % _w) + _w) % _w;
            iu1 = ((iu1 % _w) + _w) % _w;
            iv0 = ((iv0 % _h) + _h) % _h;
            iv1 = ((iv1 % _h) + _h) % _h;
        }
        //  return getValue(iu0, iv0);

        if (!_linear)
            return getValue(iu0, iv0);
        //     iu0 = 255;iu1 =256; iv0 = 599;iv1 =600;
        //   return getValue(iu0,iv0);
        return lerp(
            getValue(iu0, iv0),
            getValue(iu1, iv0),
            getValue(iu0, iv1),
            getValue(iu1, iv1),
            u,
            v);
    }

    void LoadResources() {
        LOGI("Loading {0}", _path);
        _linear = true;
        _clamp  = false;
        _isHdr  = ImageIO::isHdr(_path);
        _isRGB  = sizeof(T) == sizeof(Spectrum);
        int   w, h;
        void* pixels = nullptr;
        if (_isHdr)
            pixels = ImageIO::loadHdr(_path, TexelConversion::REQUEST_RGB, w, h).release();
        else
            pixels = ImageIO::loadLdr(_path, TexelConversion::REQUEST_RGB, w, h).release();
        if (!pixels)
            spdlog::error("Error while reading texture");
        _texels    = pixels;
        _w         = w;
        _h         = h;
        _texelType = TexelType(_isRGB << 1 | _isHdr);
        for (int y = 0; y < _h; ++y)
            for (int x = 0; x < _w; x++) {
                _average += getValue(x, y) / (float)(_w * _h);
            }
    }

    inline T getValue(int x, int y) const {
        T value;
        convertOut(_texels, value, y * _w + x, _isHdr);
        return value;
    }

    inline int width() {
        return _w;
    }
    inline int height() {
        return _h;
    }

protected:
    inline Float weight(int x, int y) const {
        if (_texelType == TexelType::RGB_HDR) {
            const vec3 v = as<vec3>()[x + y * _w];
            return max(v);
        }
        return 1;
    }

    static void convertOut(const void* data, Spectrum& value, int index, bool isHdr) {
        if (isHdr)
            value = reinterpret_cast<const vec3*>(data)[index];
        else
            value = reinterpret_cast<const Rgba*>(data)[index].normalize();
    }

    static void convertOut(const void* data, Float& value, int index, bool isHdr) {
        if (isHdr)
            value = reinterpret_cast<const Float*>(data)[index];
        else
            value = reinterpret_cast<const uint8*>(data)[index] * (1.f / 255.f);
    }

    template<typename T1>
    inline const T1* as() const {
        return reinterpret_cast<const T1*>(_texels);
    }

public:
    T average() override {
        return _average;
    }

protected:
    void*                           _texels;
    TexelType                       _texelType;
    std::unique_ptr<Distribution2D> _distribution[MAP_JACOBIAN_COUNT];
    int                             _w, _h;
    T                               _average;
    bool                            _isHdr, _isRGB, _clamp, _linear;
    std::string                     _path;
};