#include "BitMapTexture.hpp"
#include "IO/ImageIO.hpp"
#include "spdlog/spdlog.h"

static int count =0;
static vec3 rgb(0);
template<typename T>
inline const T * BitMapTexture::as() const
{
    return reinterpret_cast<const T *>(_texels);
}

inline Float BitMapTexture::weight(int x, int y) const
{   const vec3 v  =as<vec3>()[x + y*_w];
    return maxElement(v);
}

void BitMapTexture::makeSamplable(TextureMapJacobian jacobian) {
    if (_distribution[jacobian])
        return;

    std::vector<float> weights(_w*_h);
    for (int y = 0, idx = 0; y < _h; ++y) {
        float rowWeight = 1.0f;
        if (jacobian == MAP_SPHERICAL)
            rowWeight *= std::sin((y*Constant::PI)/_h);
        for (int x = 0; x < _w; ++x, ++idx)
            weights[idx] = weight(x, y)*rowWeight;
    }
//    for (int y = 0; y < _h; ++y) {
//        for (int x = 0; x < _w - 1; ++x)
//            weights[x + y*_w] = std::max(weights[x + y*_w], weights[x + 1 + y*_w]);
//        if (!_clamp)
//            weights[y*_w] = weights[_w - 1 + y*_w] = std::max(weights[_w - 1 + y*_w], weights[y*_w]);
//        for (int x = _w - 1; x > 0; --x)
//            weights[x + y*_w] = max(weights[x + y*_w], weights[x - 1 + y*_w]);
//    }
//    for (int x = 0; x < _w; ++x) {
//        for (int y = 0; y < _h - 1; ++y)
//            weights[x + y*_w] = max(weights[x + y*_w], weights[x + (y + 1)*_w]);
//        if (!_clamp)
//            weights[x] = weights[x + (_h - 1)*_w] = max(weights[x], weights[x + (_h - 1)*_w]);
//        for (int y = _h - 1; y > 0; --y)
//            weights[x + y*_w] = max(weights[x + y*_w], weights[x + (y - 1)*_w]);
//    }

    _distribution[jacobian].reset(new Distribution2D(weights.data(), _w, _h));
}

vec3 BitMapTexture::Evaluate(const Intersection * si) const {
    return vec3();
}

vec3 BitMapTexture::Evaluate(const vec2 & uv) const {
    float u = uv.x*_w;
    float v = (uv.y)*_h;
    bool linear = true;
    if (linear) {
        u -= 0.5f;
        v -= 0.5f;
    }
    int iu0 = u < 0.0f ? -int(-u) - 1 : int(u);
    int iv0 = v < 0.0f ? -int(-v) - 1 : int(v);
    int iu1 = iu0 + 1;
    int iv1 = iv0 + 1;

    u-=iu0;
    v-=iv0;
    iu0 = clamp(iu0, 0, _w - 1);
    iu1 = clamp(iu1, 0, _w - 1);
    iv0 = clamp(iv0, 0, _h - 1);
    iv1 = clamp(iv1, 0, _h - 1);

    auto ans = lerp(
            getRGB(iu0, iv0),
            getRGB(iu1, iv0),
            getRGB(iu0, iv1),
            getRGB(iu1, iv1),
            u,
            v
    );

    rgb+=ans;
    count++;

    return  ans;

}

vec2 BitMapTexture::sample(TextureMapJacobian jacobian, const vec2 & uv) const {
    vec2 newUv =_distribution[jacobian]->SampleContinuous(uv, nullptr);
    return vec2(newUv.x,1-newUv.y);
}

Float BitMapTexture::pdf(TextureMapJacobian jacobian, const vec2 & uv) const {
    vec2 newuv(uv);
    return _distribution[jacobian]->Pdf(vec2(newuv.x,(1-newuv.y))) * _w * _h;
}

vec3 BitMapTexture::getRGB(int x, int y) const {
    return as<vec3>()[x+y*_w];
}

void BitMapTexture::LoadResources( ) {
    int w, h;
    void *pixels = nullptr;
    pixels = ImageIO::loadHdr(_path, TexelConversion::REQUEST_RGB, w, h).release();
    _texels = pixels;
    _w = w;
    _h = h;
    _texelType = TexelType::RGB_HDR;

    vec3 averageRGB;
    for (int y = 0; y < _h; ++y)
        for(int x =0;x<_w;x++){
            averageRGB += getRGB(x,y)/(float)(_w * _h);
        }
    spdlog::info(toColorStr(averageRGB));
}

void BitMapTexture::debugLoginfo( ) {
    spdlog::info("\nbitMap:"+toColorStr(rgb/float(count)));
}
