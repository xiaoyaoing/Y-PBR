#include "Common/Texture.hpp"
#include "Sampler/Distrib.hpp"
#include "Colors/Spectrum.hpp"


class BitMapTexture : public  Texture<vec3>{
public:

    enum class TexelType : uint32 {
        SCALAR_LDR = 0,
        SCALAR_HDR = 1,
        RGB_LDR    = 2,
        RGB_HDR    = 3,
    };


    BitMapTexture(const std::string & path) :_path(path){

    }



    vec3 Evaluate(const Intersection * si) const override ;
    void makeSamplable(TextureMapJacobian jacobian) override ;
    vec2 sample(TextureMapJacobian jacobian, const vec2 & uv) const override;
    Float pdf(TextureMapJacobian jacobian, const vec2 & uv) const override;
    vec3 Evaluate(const vec2 & uv) const override;
    void LoadResources();
    void debugLoginfo();

protected:
    inline float weight(int x, int y) const;
    inline vec3  getRGB(int x,int y) const ;
    template<typename T>
    inline const T * as() const;

    void *_texels;
    TexelType _texelType;
    std::unique_ptr<Distribution2D> _distribution[MAP_JACOBIAN_COUNT];
    int _w,_h;
    std::string _path;
};