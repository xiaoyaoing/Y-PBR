#pragma  once

#include "Common/Texture.hpp"
#include "Sampler/Distrib.hpp"
#include "Colors/Spectrum.hpp"
#include "IO/ImageIO.hpp"
#include "stb_image.h"

struct Rgba
{
    uint8 c[4];

    vec3 normalize() const
    {
        return vec3(float(c[0]), float(c[1]), float(c[2]))*(1.0f/255.0f);
    }
};

template < class T >
class BitMapTexture : public Texture < T > {
public:

    enum class TexelType : uint32 {
        SCALAR_LDR = 0,
        SCALAR_HDR = 1,
        RGB_LDR = 2,
        RGB_HDR = 3,
    };


    BitMapTexture(const std::string & path) : _path(path) {

    }


    T Evaluate(const Intersection * si) const override {
        return Evaluate(si->uv);
    }


    void makeSamplable(TextureMapJacobian jacobian) override {
        if ( _distribution[jacobian] )
            return;

        std::vector < float > weights(_w * _h);

        for (int y = 0, idx = 0 ;y < _h ; ++ y ) {
            float rowWeight = 1.0f;
            if ( jacobian == MAP_SPHERICAL )
                rowWeight *= std::sin(( y * Constant::PI ) / _h);
            for (int x = 0 ;x < _w ;++ x, ++ idx )
                weights[idx] =weight(x, y) * rowWeight;
        }
        _distribution[jacobian].reset(new Distribution2D(weights.data(), _w, _h));
    }

    vec2 sample(TextureMapJacobian jacobian, const vec2 & uv, Float * pdf) const override {
        vec2 newUv = _distribution[jacobian]->SampleContinuous(uv, pdf);
        return vec2(newUv.x, newUv.y);
    }

    Float pdf(TextureMapJacobian jacobian, const vec2 & uv) const override {
        vec2 newuv(uv);
        Float res = _distribution[jacobian]->Pdf(vec2(newuv.x,  newuv.y ));
        if( isnan(res) || res == 0){

        }
        return res;
    }

    T Evaluate(const vec2 & uv) const override {
        float u = uv.x * _w;
        float v = uv.y * _h;
        bool linear = true;
        if ( linear ) {
            u -= 0.5f;
            v -= 0.5f;
        }
        int iu0 = u < 0.0f ? - int(- u) - 1 : int(u);
        int iv0 = v < 0.0f ? - int(- v) - 1 : int(v);
        int iu1 = iu0 + 1;
        int iv1 = iv0 + 1;

        u -=iu0;
        v -=iv0;
        iu0 = clamp(iu0, 0, _w - 1);
        iu1 = clamp(iu1, 0, _w - 1);
        iv0 = clamp(iv0, 0, _h - 1);
        iv1 = clamp(iv1, 0, _h - 1);

        auto ans = lerp(
                getValue(iu0, iv0),
                getValue(iu1, iv0),
                getValue(iu0, iv1),
                getValue(iu1, iv1),
                u,
                v
        );
        return ans;
    }

    void LoadResources(){
        isHdr = stbi_is_hdr(_path.c_str());
        isRgb = sizeof(T) == sizeof(Spectrum);
        int w, h;
        void * pixels = nullptr;
        if ( isHdr )
            pixels = ImageIO::loadHdr(_path, TexelConversion::REQUEST_RGB, w, h).release();
        else
            pixels = ImageIO::loadLdr(_path, TexelConversion::REQUEST_RGB, w, h).release();
        _texels = pixels;
        _w = w;
        _h = h;
        _texelType = TexelType(isRgb<<1 | isHdr);
        for (int y = 0 ; y < _h ; ++ y )
            for (int x = 0 ;x < _w ;x ++ ) {
                _average +=getValue(x, y) / (float) (_w * _h);
        }
    }

protected:

    inline Float weight(int x, int y) const {
        if ( _texelType == TexelType::RGB_HDR ) {
            const vec3 v = as < vec3 >()[x + y * _w];
            return luminace(v);
        }
    }

    inline T getValue(int x, int y) const {
        T value;
        convertOut(_texels,value,y*_w + x,isHdr);
        return value;
    }

    static void convertOut(const void * data ,Spectrum & value,int index,bool isHdr){
        if(isHdr)
            value = reinterpret_cast<const vec3 *>(data)[index];
        else
            value = reinterpret_cast<const Rgba *>(data)[index].normalize();
    }

    static void convertOut(const void * data ,Float & value,int index,bool isHdr){
        if(isHdr)
            value = reinterpret_cast<const Float *>(data)[index];
        else
            value =  reinterpret_cast<const uint8 *>(data)[index] * (1.f/255.f);
    }

    template < typename T1 >
    inline const T1 * as( ) const {
        return reinterpret_cast<const T1 *>(_texels);
    }

public:
    T average( ) override {
        return _average;
    }

protected:
    void * _texels;
    TexelType _texelType;
    std::unique_ptr < Distribution2D > _distribution[MAP_JACOBIAN_COUNT];
    int _w, _h;
    T _average;
    bool isHdr;
    bool isRgb;
    std::string _path;
};