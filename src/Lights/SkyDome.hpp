#include "Light.hpp"
#include "Common/Json.hpp"
#include "Texture/BitMapTexture.hpp"
class SkyDome : public  Light{
public:
    SkyDome(const Json & json);

    Spectrum
    sampleLi(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf, VisibilityTester * vis) const override;

    LightSampleResult sampleDirect(const vec2 & positionSample, const vec2 & dirSample) override;


    Spectrum Le(const Ray & ray) const override;

    Spectrum Power( ) override;

    void Preprocess(const Scene & scene) override;

    Float PdfLi(const Intersection & pShape, const vec3 & ref) const override;

    bool isDeltaLight( ) const override;

protected:
    vec2 directionToUV(const vec3 &wi) const;
    vec2 directionToUV(const vec3 &wi, float &sinTheta) const;
    vec3 uvToDirection(vec2 uv, float &sinTheta) const;

    std::shared_ptr<BitMapTexture<Spectrum>> _sky;
    float _temperature;
    float _gammaScale;
    float _turbidity;
    float _intensity;
    bool _doSample;
    mat4  _transform;

    vec3 _worldCenter;
    Float _worldRadius;
};