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
    return luminace(v);
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

    _distribution[jacobian].reset(new Distribution2D(weights.data(), _w, _h));
}

vec3 BitMapTexture::Evaluate(const Intersection * si) const {
    return vec3();
}

vec3 BitMapTexture::Evaluate(const vec2 & uv) const {
    float u = uv.x *_w;
    float v = uv.y *_h;
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
    return  ans;

}

vec2 BitMapTexture::sample(TextureMapJacobian jacobian, const vec2 & uv,Float * pdf) const {
    vec2 newUv =_distribution[jacobian]->SampleContinuous(uv,  pdf);
    return vec2(newUv.x,newUv.y);
}

Float BitMapTexture::pdf(TextureMapJacobian jacobian, const vec2 & uv) const {
    vec2 newuv(uv);
    Float res =  _distribution[jacobian]->Pdf(vec2(newuv.x,(newuv.y)));
    return res;
    return _distribution[jacobian]->Pdf(vec2(newuv.x,(newuv.y))) * _w * _h;
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

    for (int y = 0; y < _h; ++y)
        for(int x =0;x<_w;x++){
           _average += getRGB(x,y)/(float)(_w * _h);
    }
}

void BitMapTexture::debugLoginfo( ) {
    spdlog::info("\nbitMap:"+toColorStr(rgb/float(count)));
}

vec3 BitMapTexture::average( ) {
    return _average;
}
